# Collatz Conjecture tester

This optimized [Collatz Conjecture](https://en.wikipedia.org/wiki/Collatz_conjecture) tester can test the provided numbers for how many iterations they need before they reach the number 1. The path is not printed, only statistics: the number of 3*n+1 ("mul3") and n/2 ("div2") operations, and total. As a progress indicator prints the current length of the number in bytes (every 25000 iterations).

If you're lucky to find a number that proves that Collatz Conjecture is false (or a bug in the program), then the program will hang in an infinite loop.

* Note: Doesn't work with negative numbers.

## Building

1. `make all`  
default build
2. `make profile PROF_ARGS="--ones 1000000"`  
build with the profile-guided optimizations using GCC (add CLANG=1 when using Clang)

* `make check` for a simple build check.
* Compile-time option `NBITS`, greatly affects performance, should be the number of bits in the machine register. The default is 64 for 64-bit machines (also needs the __int64 type supported by the compiler), otherwise 32. 

## Usage

`collatz_test [--lut 1..26] [mode] {num|file}`

## Options

`--lut N` builds a lookup table to speed up the process by a factor of N (default 20). LUT consumes `2^n * 3 * NBITS/8` bytes of memory. The LUT depth is also limited to `5 * NBITS/8`.

* Doesn't affect performance in the latest version with dynamic expansion of LUT values.

## Modes

`--num` read decimal number from command line (default)  
`--file` load specified file as a little endian number  
`--ones` test 2^n-1  

## Example

	$ ./collatz_test 989345275647
	mul3 = 506, div2 = 842, total = 1348
	$ ./collatz_test 0xffffffff
	mul3 = 162, div2 = 289, total = 451
	$ ./collatz_test --ones 1000000
	lut: 0.016s
	bytes: 198128
	bytes: 145760
	bytes: 94000
	bytes: 41896
	mul3 = 4805005, div2 = 8615753, total = 13420758
	time: 1.155s
	$ head -c100000 /dev/urandom > random.bin
	$ ./collatz_test --file random.bin
	lut: 0.016s
	bytes: 47936
	mul3 = 1923602, div2 = 3848837, total = 5772439
	time: 0.213s

