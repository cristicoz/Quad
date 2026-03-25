Quad (Sabs.Numerics.Quad)
Copyright (c) 2024 Cristian Cozma. Licensed under MIT License.
A high-performance, lightweight library for Double-Double arithmetic, providing 106-bit precision within the .NET ecosystem.
Unlike traditional implementations, Quad introduces a specialized class of numerical kernels: Fused Error-Free Transformations (FEFT). This approach achieves maximum precision using a minimal number of atomic operations (e.g., Double-Double addition in just 10 operations), specifically optimized for the execution pipelines of modern CPUs.
Why Quad?
Native Speed: Significantly faster than System.Decimal by operating directly on FPU/SSE registers.
106-Bit Precision: Manages ~32 significant digits, capturing the full precision of any built-in C# data type, including long and decimal.
Financial Range: Can sum massive amounts of decimal data without the overflow limit, utilizing the double exponent (
).
Robustness: Maintains accuracy even in the denormalized range down to Epsilon (
).
Technical Validation (Stress Tests)
Test Case	Input	Round-trip Result	Status
Max Value	Log2(MaxValue) -> Pow2()	Bit-identical (h/l)	PASSED
Decimal Max	(decimal)(Quad)decimal.MaxValue	79228162514264337593543950335	PASSED
Subnormal	Pow2(Log2(Epsilon * 3))	Bit-identical (h/l)	PASSED
Boundary	long.MaxValue + 0.999...16	Preserves 106-bit limit	PASSED
The FEFT Advantage
At its core, Quad replaces "cosmetic" normalization with raw error-bit conservation. The 10-operation DWPlusDW kernel creates an aggressive self-correcting accumulator. This makes it the ideal engine for N-body physical simulations, massive financial aggregations, and high-fidelity engineering where performance and bit-integrity are paramount.
