# Instruction for measuring the performance of Hashing

On terminal 1
```
nohup ./server_alpha 2147480000 1677762560 10240 32 100000 > Hashing_log_26_08_21_06.txt 2>&1 
```

On terminal 2
```
./server_beta
```

On terminal 3
```
./server_gamma
```