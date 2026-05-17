# Bamboo Filter

C++ implementation of the Bamboo Filter (Wang et al. 2022) — a resizable approximate membership query filter that extends Cuckoo Filters with smooth, incremental expansion and shrinkage. Built for the FER course *Bioinformatika 1* (2025./2026.), task (2).

**Goal.** Search for random k-mers (k ∈ {10, 20, 50, 100, 200}) in the *E. coli* K-12 MG1655 genome and in synthetic DNA, using our own Bamboo Filter, and compare against the original reference implementation.

**Authors:** Antonio Šimić, Jakov Malić

## References

- Wang et al. *Bamboo Filters: Make Resizing Smooth* (ICDE 2022)
- Wang et al. *Bamboo Filters: Make Resizing Smooth and Adaptive* (IEEE TNET 2024)
- Fan et al. *Cuckoo Filter: Practically Better Than Bloom* (CoNEXT 2014)
- Reference implementation: <https://github.com/wanghanchengchn/bamboofilters>
