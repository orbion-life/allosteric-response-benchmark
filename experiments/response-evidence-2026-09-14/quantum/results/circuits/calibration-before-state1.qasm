OPENQASM 3.0;
bit[6] b;
#pragma braket verbatim
box {
  prx(3.1415926535897931, 0) $7;
  prx(3.1415926535897931, 0) $8;
  prx(3.1415926535897931, 0) $9;
  prx(3.1415926535897931, 0) $10;
  prx(3.1415926535897931, 0) $13;
  prx(3.1415926535897931, 0) $14;
}
b[0] = measure $7;
b[1] = measure $8;
b[2] = measure $9;
b[3] = measure $10;
b[4] = measure $13;
b[5] = measure $14;
