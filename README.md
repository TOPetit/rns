# RNS
Internship repository

# Config
```shell
sudo modprobe msr
sudo /bin/echo "-1" > /proc/sys/kernel/perf_event_paranoid
sudo echo 2 | dd of=/sys/devices/cpu/rdpmc
sudo echo 2 | dd of=/sys/bus/event_source/devices/cpu/rdpmc
sudo wrmsr -a 0x38d 0x0333
```
# TODO

## Timings:
1. Instructions per cycles

## Vectorisation:
1. Base conversion
2. Modular multplication

## Fix:
1. Vectorized substraction