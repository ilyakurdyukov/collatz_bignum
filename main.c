/*
 * Copyright (C) 2021 Ilya Kurdyukov
 *
 * Optimized Collatz Conjecture tester.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#ifdef __SSE2__
#include <emmintrin.h>
#endif

#ifdef __SSSE3__
#include <tmmintrin.h>
#else
#define _mm_alignr_epi8(a, b, n) \
	_mm_or_si128(_mm_bsrli_si128(b, n), _mm_bslli_si128(a, 16 - n))
#endif

#ifdef __AVX2__
#include <immintrin.h>
#endif

#if 1 && defined(__GNUC__)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)
#define LIKELY(x) __builtin_expect(!!(x), 1)
#else
#define UNLIKELY(x) x
#define LIKELY(x) x
#endif

#include <time.h>
#include <sys/time.h>
static int64_t get_time_usec() {
	struct timeval time;
	gettimeofday(&time, NULL);
	return time.tv_sec * (int64_t)1000000 + time.tv_usec;
}

#ifdef NBITS
#define N NBITS
#else
#ifdef __LP64__
#define N 64
#else
#define N 32
#endif
#endif

#if N == 64
#define T uint64_t
#define T2 unsigned __int128
#elif N == 32
#define T uint32_t
#define T2 uint64_t
#elif N == 16
#define T uint16_t
#define T2 uint32_t
#elif N == 8
#define T uint8_t
#define T2 uint16_t
#else
#error
#endif

typedef struct {
	size_t cur, max, inc; T *buf;
} bignum_t;

#define BLOCKBYTES 4096

static void bignum_init(bignum_t *bn) {
	bn->cur = 0;
	bn->max = bn->inc = BLOCKBYTES / sizeof(T);
	bn->buf = (T*)malloc(BLOCKBYTES);
}

static void* bignum_alloc(bignum_t *bn, size_t n) {
	size_t align = BLOCKBYTES, max = (n + align - 1) & -align;
	uint8_t *buf = (uint8_t*)malloc(max);
	bn->inc = BLOCKBYTES / sizeof(T);
	bn->cur = 0;
	bn->max = (max + sizeof(T) - 1) / sizeof(T);
	bn->buf = (T*)buf;
	return buf;
}

static void bignum_prepare(bignum_t *bn, size_t n) {
	union { T i; uint8_t b; } check = { 1 };
	size_t i, k;
	uint8_t *buf = (uint8_t*)bn->buf;
	while (n > 1 && !buf[n - 1]) n--;
	k = n & (sizeof(T) - 1);
	if (k) memset(buf + n, 0, sizeof(T) - k);
	bn->cur = (n + sizeof(T) - 1) / sizeof(T);
	// swap for big-endian CPUs
	if (!check.b)
	for (k = 0; k < n; k += sizeof(T))
	for (i = 0; i < sizeof(T) / 2; i++) {
		uint8_t a = buf[k + i], b = buf[k + sizeof(T) - 1 - i];
		buf[k + i] = b;
		buf[k + sizeof(T) - 1 - i] = a;
	}
}

static void* bignum_read(bignum_t *bn, const char *fn) {
	size_t n, j;
	FILE *fi = fopen(fn, "rb");
	bn->buf = NULL;
	if (fi) {
		fseek(fi, 0, SEEK_END);
		n = ftell(fi);
		fseek(fi, 0, SEEK_SET);
		if (n) {
			uint8_t *buf = (uint8_t*)bignum_alloc(bn, n);
			if (buf) {
				j = fread(buf, 1, n, fi);
				bignum_prepare(bn, j);
			}
		}
		fclose(fi);
	}
	return bn->buf;
}

#define bignum_mul_add(x, a, b) bignum_shrN_mul_add(x, 0, a, b)
static void bignum_shrN_mul_add(bignum_t *bn, size_t shr, T mul, T add) {
	size_t i = shr / N, j = 0, n = bn->cur;
	T *buf = bn->buf;
	unsigned k = shr & (N - 1);
	if (LIKELY(i < n)) {
#if 1 && defined(__AVX2__) && N == 32
#ifdef __GNUC__
// _mm256_zextsi128_si256 doesn't supported in GCC < 10
// in this code v3 can't have a dirty high half
#if __GNUC__ < 10
#define _mm256_zextsi128_si256 _mm256_castsi128_si256
#endif
#endif
		if (LIKELY(i + 7 < n)) {
#define O 1ULL << 63
#define I O | 1
			static const uint64_t __attribute__((aligned(32))) carry_tab[16][4] = {
				{O,O,O,O}, {I,O,O,O}, {O,I,O,O}, {I,I,O,O},
				{O,O,I,O}, {I,O,I,O}, {O,I,I,O}, {I,I,I,O},
				{O,O,O,I}, {I,O,O,I}, {O,I,O,I}, {I,I,O,I},
				{O,O,I,I}, {I,O,I,I}, {O,I,I,I}, {I,I,I,I} };
#undef I
#undef O
			__m128i v3 = _mm_cvtsi32_si128(add);
			__m256i v0, v1, v2 = _mm256_set1_epi32(mul), v4, v5;
			__m256i c1 = _mm256_set1_epi32(-1), c2 = _mm256_slli_epi64(c1, 63);
			int max_mask, carry = 0;
			v3 = _mm_sll_epi64(v3, _mm_cvtsi32_si128(k));
			c1 = _mm256_srli_epi64(c1, 1);

			v4 = _mm256_loadu_si256((__m256i*)&buf[i]); i += 8;
			v0 = _mm256_zextsi128_si256(_mm_cvtsi32_si128((1 << k) - 1));
			v1 = _mm256_mul_epu32(_mm256_bsrli_epi128(v4, 4), v2);
			v4 = _mm256_andnot_si256(v0, v4);
			v0 = _mm256_add_epi64(_mm256_zextsi128_si256(v3), _mm256_mul_epu32(v4, v2));

#define PART1(p1, p2) \
	v4 = _mm256_permute2x128_si256(v1, p1, p2); v5 = v1; \
	v1 = _mm256_alignr_epi8(v1, v4, 12); \
	v0 = _mm256_xor_si256(v0, c2); \
	v1 = _mm256_add_epi64(v1, v0); \
	v0 = _mm256_cmpgt_epi64(v0, v1); \
	v4 = _mm256_cmpeq_epi64(v1, c1);
#define PART2 \
	max_mask = _mm256_movemask_pd(_mm256_castsi256_pd(v4)); \
	carry += _mm256_movemask_pd(_mm256_castsi256_pd(v0)) * 2; \
	carry += max_mask; max_mask ^= carry; carry >>= 4; \
	v0 = _mm256_load_si256((__m256i*)carry_tab[max_mask & 15]); \
	v0 = _mm256_add_epi32(v1, v0);

			PART1(v1, 0x08)
			while (i + 7 < n) {
				PART2
				v1 = _mm256_loadu_si256((__m256i*)&buf[i]); i += 8;
				_mm256_storeu_si256((__m256i*)&buf[j], v0); j += 8;
				v0 = _mm256_mul_epu32(v1, v2);
				v1 = _mm256_mul_epu32(_mm256_bsrli_epi128(v1, 4), v2);
				PART1(v5, 0x03)
			}
			PART2
			_mm256_storeu_si256((__m256i*)&buf[j], v0); j += 8;
			v3 = _mm256_extracti128_si256(v5, 1);
			add = _mm_extract_epi32(v3, 3) + carry;
#undef PART2
#undef PART1
		} else
#elif 1 && defined(__SSE2__) && N == 32
		if (LIKELY(i + 3 < n)) {
#define O 1U << 31
#define I O | 1
			static const uint32_t __attribute__((aligned(16))) carry_tab[16][4] = {
				{O,O,O,O}, {I,O,O,O}, {O,I,O,O}, {I,I,O,O},
				{O,O,I,O}, {I,O,I,O}, {O,I,I,O}, {I,I,I,O},
				{O,O,O,I}, {I,O,O,I}, {O,I,O,I}, {I,I,O,I},
				{O,O,I,I}, {I,O,I,I}, {O,I,I,I}, {I,I,I,I} };
#undef I
#undef O
			__m128i v3 = _mm_cvtsi32_si128(add);
			__m128i v0, v1, v2 = _mm_set1_epi32(mul), v4, v5;
			__m128i c1 = _mm_set1_epi32(-1), c2 = _mm_slli_epi32(c1, 31);
			int max_mask, carry = 0;
			v3 = _mm_sll_epi64(v3, _mm_cvtsi32_si128(k));
			c1 = _mm_srli_epi32(c1, 1);

			v4 = _mm_loadu_si128((__m128i*)&buf[i]); i += 4;
			v0 = _mm_cvtsi32_si128((1 << k) - 1);
			v1 = _mm_mul_epu32(_mm_bsrli_si128(v4, 4), v2);
			v4 = _mm_andnot_si128(v0, v4);
			v0 = _mm_add_epi64(v3, _mm_mul_epu32(v4, v2));

#define PART1(p1) \
	v4 = p1; v5 = v1; \
	v0 = _mm_xor_si128(v0, c2); \
	v1 = _mm_add_epi32(v4, v0); \
	v0 = _mm_cmpgt_epi32(v0, v1); \
	v4 = _mm_cmpeq_epi32(v1, c1);
#define PART2 \
	max_mask = _mm_movemask_ps(_mm_castsi128_ps(v4)); \
	carry += _mm_movemask_ps(_mm_castsi128_ps(v0)) * 2; \
	carry += max_mask; max_mask ^= carry; carry >>= 4; \
	v0 = _mm_load_si128((__m128i*)carry_tab[max_mask & 15]); \
	v0 = _mm_add_epi32(v1, v0);

			PART1(_mm_bslli_si128(v1, 4))
			while (i + 3 < n) {
				PART2
				v1 = _mm_loadu_si128((__m128i*)&buf[i]); i += 4;
				_mm_storeu_si128((__m128i*)&buf[j], v0); j += 4;
				v0 = _mm_mul_epu32(v1, v2);
				v1 = _mm_mul_epu32(_mm_bsrli_si128(v1, 4), v2);
				PART1(_mm_alignr_epi8(v1, v5, 12))
			}
			PART2
			_mm_storeu_si128((__m128i*)&buf[j], v0); j += 4;
			add = _mm_cvtsi128_si32(_mm_bsrli_si128(v5, 12)) + carry;
#undef PART2
#undef PART1
		} else
#endif
		{
			// note: shift optimization leaves (shr % N) binary zeros at the beginning
			T2 tmp = (buf[i++] & (T)((T)-1 << k)) * (T2)mul + ((T2)add << k);
			buf[j++] = (T)tmp;
			add = tmp >> N;
		}
		while (LIKELY(i < n)) {
			T2 tmp = buf[i++] * (T2)mul + add;
			buf[j++] = (T)tmp;
			add = tmp >> N;
		}
	}
	if (LIKELY(add)) {
		size_t max = bn->max;
		if (UNLIKELY(j >= max)) {
			bn->max = max += bn->inc;
			bn->buf = buf = (T*)realloc(buf, max * sizeof(T));
			if (UNLIKELY(!buf)) {
				printf("!!! realloc failed\n");
				exit(1);
			}
		}
		buf[j++] = add;
	}
	bn->cur = j;
}

static size_t bignum_ctz(bignum_t *bn, int *end) {
	size_t i = 0, k = 0, n = bn->cur;
	T a = 0, *buf = bn->buf;
	while (LIKELY(i < n)) {
		a = buf[i++];
		if (LIKELY(a)) {
#if 1 && defined(__GNUC__)
			T x = i < n ? a & ~(a >> 2 & ~(a >> 1)) : a; /* 101 ^ 001 */
			a >>= k = N <= sizeof(int) ? __builtin_ctz(x) : __builtin_ctzll(x);
			k = N - k;
#else
			for (k = N; !(a & 1); k--) a >>= 1;
			// shortcut: reduce trailing 1(01) sequences to 1
			// so this function is no longer a simple CTZ
			if (i < n) while ((a & 7) == 5) { a >>= 2; k -= 2; }
#endif
			break;
		}
	}
	*end = !(a >> 1 | (i - n));
	return i * N - k;
}

int main(int argc, char **argv) {
	bignum_t bn; int lut = 20;
	typedef struct { T mul, add, inc; } lut_t;
	static const lut_t lut1 = { 3, 2, 1 };
	lut_t *lut_ptr = (lut_t*)&lut1;

	if (argc > 2 && !strcmp(argv[1], "--lut")) {
		lut = atoi(argv[2]);
		if (lut < 1) lut = 1;
		if (lut > 26) lut = 26;
		argc -= 2; argv += 2;
	}

	if ((argc | 1) != 3) {
		printf(
			"Usage:\n"
			"  collatz_test [--lut 1..26] [mode] {num|file}\n"
			"Modes:\n"
			"  --num   read decimal/hex number from command line (default)\n"
			"  --file  load specified file as a little endian number\n"
			"  --ones  test 2^n-1\n"
		);
		return 1;
	}

	if (argc == 2 || !strcmp(argv[1], "--num")) {
		const char *s = argv[argc - 1];
		unsigned a, c, base = 10;
		bignum_init(&bn);
		*bn.buf = 0;
		if (*s == '0' && (s[1] == 'X' || s[1] == 'x')) {
			base = 16; s += 2;
		}
		while ((c = *s++)) {
			if ((a = c - '0') > 9) {
#if 0 /* portable */
				const char *hex = "AaBbCcDdEeFf", *ptr;
				if (base <= 10 || !(ptr = strchr(hex, c))) {
					printf("!!! unexpected character in a number string\n");
					return 1;
				}
				a = 10 + ((ptr - hex) >> 1);
#else /* ASCII specific */
				if (base <= 10 || (a = (c | 0x20) - 'a') > 5) {
					printf("!!! unexpected character in a number string\n");
					return 1;
				}
				a += 10;
#endif
			}
			bignum_mul_add(&bn, base, a);
		}
	} else if (!strcmp(argv[1], "--file")) {
		if (!bignum_read(&bn, argv[2])) {
			printf("!!! bignum_read failed\n");
			return 1;
		}
	} else if (!strcmp(argv[1], "--ones")) {
		int n = atoi(argv[2]), bytes = (n + 7) >> 3;
		if (n < 1) return 1;
		if (!bignum_alloc(&bn, bytes)) {
			printf("!!! bignum_alloc failed\n");
			return 1;
		}
		memset(bn.buf, -1, bytes);
		if (n & 7) ((uint8_t*)bn.buf)[bytes - 1] = (1 << (n & 7)) - 1;
		bignum_prepare(&bn, bytes);
	} else {
		printf("!!! unknown mode\n");
		return 1;
	}

	if (bn.cur < 3) lut = 1;
	if (lut > (N / 8) * 5) lut = (N / 8) * 5;

	if (lut > 1) {
		int64_t t0 = get_time_usec(), t1;
		int i, j, n = 1 << lut;
		lut_ptr = malloc((n >> 1) * sizeof(lut_t));
		for (i = 1; i < n; i += 2) {
			T a = i, m = 1; int n = 0;
			for (j = 0; j < lut; j++) {
#if 0
				a = a & 1 ? n++, m *= 3, (a >> 1) + a + 1 : a >> 1;
#else /* branchless */
				T x = a & 1;
				n += x; m += (m << 1) & -x;
				a = (a >> 1) + ((a + 1) & -x);
#endif
			}
			lut_ptr[i >> 1].mul = m;
			lut_ptr[i >> 1].add = a;
			lut_ptr[i >> 1].inc = n;
		}
		t1 = get_time_usec();
		if (t1 - t0 > 10000)
			printf("lut: %.3fs\n", (t1 - t0) * 0.000001);
	}

	{
		long long mul3 = 0, div2 = 0;
		int64_t t0 = get_time_usec(), t1;
		int i = 0;

		for (;;) {
			int end; size_t shr, pos;
			int inc = 1; T mul = 3, add = 1;

			div2 += shr = bignum_ctz(&bn, &end);
			pos = shr / N;
			if (UNLIKELY(end)) break;
			if (LIKELY(pos + 2 < bn.cur)) {
				unsigned k = shr & (N - 1); size_t i;
				T bits = bn.buf[pos++] >> k;
				bits |= (T2)bn.buf[pos] << (N - k);
				i = (bits & ((1 << lut) - 1)) >> 1;
				mul = lut_ptr[i].mul;
				add = lut_ptr[i].add;
				inc = lut_ptr[i].inc;
#if 0
				div2 += lut; shr += lut;
#else
				// extend values from LUT
				bits >>= i = lut;
				// same condition as: mul <= (T)-1 / 3
				while (LIKELY(inc < (N / 8) * 5)) {
					T x;
					add += mul & -(bits & 1); bits >>= 1;
					x = add & 1;
					inc += x; mul += (mul << 1) & -x;
					add = (add >> 1) + ((add + 1) & -x);
					if (UNLIKELY(!(++i & (N - 1)))) {
						if (UNLIKELY(pos + 2 >= bn.cur)) break;
						bits = bn.buf[pos++] >> k;
						bits |= (T2)bn.buf[pos] << (N - k);
					}
				}
				div2 += i; shr += i;
#endif
			}
			mul3 += inc;
			div2 -= shr & (N - 1);	// compensation
			bignum_shrN_mul_add(&bn, shr, mul, add);
			if (UNLIKELY(++i >= 25000)) {
				i = 0;
				printf("bytes: %ld\n", (long)bn.cur * sizeof(T));
			}
		}
		printf("mul3 = %lld, div2 = %lld, total = %lld\n", mul3, div2, mul3 + div2);
		t1 = get_time_usec();
		if (t1 - t0 > 10000)
			printf("time: %.3fs\n", (t1 - t0) * 0.000001);
	}
	return 0;
}
