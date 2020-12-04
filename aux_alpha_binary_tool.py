#!/usr/bin/python3

alphanumerics = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ+-*/%"
assert(len(alphanumerics) == 41)

def int_to_alphanumeric(x:int, l:int=None) -> str:
    digits = []
    while True:
        digits.append(x % len(alphanumerics))
        x = x // len(alphanumerics)
        if x == 0: break
    
    x_str = "".join([alphanumerics[i] for i in reversed(digits)])
    if l:
        x_str = x_str.rjust(l, alphanumerics[0])
    return x_str

def alphanumeric_to_int(a:str) -> int:
    a = reversed(a)
    res = 0
    factor = 1
    for c in a:
        res += alphanumerics.index(c)*factor
        factor *= len(alphanumerics)
    return res

def create_table_2byte_to_3alphanumeric(file_name:str=None):
    out_str = ""
    for i in range(256**2):
        out_str += f'"{int_to_alphanumeric(i, 3)}", '

    if file_name:
        out_file = open(file_name, "w+")
        out_file.write(out_str)
        out_file.close()
    else:
        print(out_str)

def create_table_3alphanumeric_to_2bytes(file_name:str=None):
    out_dict = {}
    for a in alphanumerics:
        for b in alphanumerics:
            for c in alphanumerics:
                alpha = f'{a}{b}{c}'
                out_dict[alpha] = alphanumeric_to_int(alpha)


    if file_name:
        out_file = open(file_name, "w+")
        out_file.write(out_str)
        out_file.close()
    else:
        print(out_str)

def test_conversions():
    for i in range(256**2):
        alpha = int_to_alphanumeric(i, 3)
        res_i = alphanumeric_to_int(alpha)
        if not i == res_i:
            print(f'Wrong {i} -> {alpha} -> {res_i}')
            return

    for a in alphanumerics:
        for b in alphanumerics:
            for c in alphanumerics:
                alpha = f'{a}{b}{c}'
                i = alphanumeric_to_int(alpha)
                res_alpha = int_to_alphanumeric(i, 3)
                if not alpha == res_alpha:
                    print(f'Wrong {alpha} -> {i} -> {res_alpha}')
    print("Done")

out_str = ""
for c in alphanumerics:
    out_str += f'\'{c}\','
print(out_str)