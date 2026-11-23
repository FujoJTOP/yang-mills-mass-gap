import struct, sys

base, out, module = sys.argv[1], sys.argv[2], sys.argv[3]
L, N, NLINK = int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
betas = [float(x) for x in sys.argv[7].split(',')]

src = open(base).read()
src = src.replace('module z2b', f'module {module}')
src = src.replace('const L: u32 = 4;', f'const L: u32 = {L};')
src = src.replace('const N: u32 = 256;', f'const N: u32 = {N};')
src = src.replace('const NLINK: u32 = 1024;', f'const NLINK: u32 = {NLINK};')

s = src.index('fn beta_bits')
e = src.index('fn main')
new = 'fn beta_bits(c: u32) -> u64 {\n'
for i, x in enumerate(betas):
    b = hex(struct.unpack('>Q', struct.pack('>d', x))[0])
    if i < len(betas) - 1:
        new += f'    if c == {i} {{ return {b}; }}\n'
    else:
        new += f'    return {b};\n'
new += '}\n\n'
src = src[:s] + new + src[e:]
open(out, 'w').write(src)
print(out, 'L', L, 'N', N, 'NLINK', NLINK)
print('betas', betas)
