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

#define bignum_mul_add(x, a, b) bignum_shr_mul_add(x, 0, a, b)
static void bignum_shr_mul_add(bignum_t *bn, size_t shr, T mul, T add) {
	size_t i = shr / N, j = 0, n = bn->cur;
	T *buf = bn->buf, x;
	unsigned k = shr & (N - 1);
	x = buf[i++] >> k;
	while (LIKELY(i < n)) {
		T2 tmp = (T2)buf[i++] << (N - k) | x;
		x = tmp >> N;
		tmp = (T)tmp * (T2)mul + add;
		buf[j++] = (T)tmp;
		add = tmp >> N;
	}
	if (LIKELY(x)) {
		T2 tmp = x * (T2)mul + add;
		buf[j++] = (T)tmp;
		add = tmp >> N;
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
			a >>= k = N <= sizeof(int) ? __builtin_ctz(a) : __builtin_ctzll(a);
			k = N - k;
#else
			for (k = N; !(a & 1); k--) a >>= 1;
#endif
			// shortcut: reduce trailing 1(01) sequences to 1
			// so this function is no longer a simple CTZ
			if (i < n) while ((a & 7) == 5) { a >>= 2; k -= 2; }
			break;
		}
	}
	*end = a < 2 && i == n;
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

	if (argc < 3) {
		printf(
			"Usage:\n"
			"  collatz_test [--lut 1..26] mode {num|file}\n"
			"Modes:\n"
			"  --num   read decimal number from command line\n"
			"  --file  load specified file as a little endian number\n"
			"  --ones  test 2^n-1\n"
		);
		return 1;
	}

	if (!strcmp(argv[1], "--num")) {
		const char *s = argv[2];
		bignum_init(&bn);
		*bn.buf = 0;
		while (*s) {
			char a = *s++ - '0';
			if ((unsigned)a > 9) {
				printf("!!! unexpected character in a decimal number string\n");
				return 1;
			}
			bignum_mul_add(&bn, 10, a);
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
			T x; int end; size_t shr;
			int inc = 1; T mul = 3, add = 1;

			div2 += shr = bignum_ctz(&bn, &end);
			if (UNLIKELY(end)) break;
			{
				size_t i = shr / N;
				unsigned k = shr & (N - 1);
				x = bn.buf[i++] >> k;
				if (LIKELY(i < bn.cur)) x |= (T2)bn.buf[i] << (N - k);
			}
			x = (x & ((1 << lut) - 1)) >> 1;

			if (LIKELY(shr / N + 2 < bn.cur)) {
				div2 += lut; shr += lut;
				mul = lut_ptr[x].mul;
				add = lut_ptr[x].add;
				inc = lut_ptr[x].inc;
			}
			mul3 += inc;
			bignum_shr_mul_add(&bn, shr, mul, add);
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
