# RNS
2021 internship repository.
> A residue numeral system (RNS) is a numeral system representing integers by their values modulo several pairwise coprime integers called the moduli. This representation is allowed by the Chinese remainder theorem, which asserts that, if N is the product of the moduli, there is, in an interval of length N, exactly one integer having any given set of modular values. The arithmetic of a residue numeral system is also called multi-modular arithmetic.

Wikipedia.

# Config
```shell
sudo modprobe msr
sudo /bin/echo "-1" > /proc/sys/kernel/perf_event_paranoid
sudo echo 2 | dd of=/sys/devices/cpu/rdpmc
sudo echo 2 | dd of=/sys/bus/event_source/devices/cpu/rdpmc
sudo wrmsr -a 0x38d 0x0333
```
or
```shell
chmod +x ./init.sh
sudo ./init.sh
```
# TODO

## Timings:
- [x] Instructions per cycles
- [x] Explain cycles for first call

## Vectorisation:
- [x] Base conversion
- [x] Modular multplication
- [x] Cox modular multiplication
- [x] Cox base extension
- [x] Cox computing k
- [ ] ~~Hierarchical base extension~~

## Fix:
- [x] Vectorized substraction
- [x] Better m256i -> int64_t rns conversion

## Test:
- [ ] Conversion
- [ ] Addition
- [ ] Multiplication
- [ ] Substraction
- [ ] Base conversion
- [ ] Modular multiplication