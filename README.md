# fountain
Fountain Codes for QR transmission

```
make
./main <chunk byte length> <num base chunks>
```

'base chunks' are the real data to transfer, each of same byte length.
Example:

```
./main 1024 50
```
For 50KB total data

