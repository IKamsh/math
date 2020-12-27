def reduce_bin_form(a, b, c):
    '''
    reduce quadratic positive difinite binary form
    ax^2 + bxy + cy^2
    (positive difinite means that a and c >= 0)
    '''
    
    def clip_b1(a, b):
        assert abs(b) > a
        b1 = lambda a,b,u: b + 2 * u * a
        if b < 0:
            u = 1
            while abs(b1(a, b, u)) > a:
                u += 1
        elif b > 0:
            u = -1
            while abs(b1(a,b,u)) > a:
                u -= 1
        return b1(a,b,u)
    
    if c < a:
        return reduce_bin_form(c, -b, a)
    elif abs(b) > a:
        b1 = clip_b1(a, b)
        d = b * b - 4 * a * c
        c1 = (b1 * b1 - d)//(4 * a)
        return reduce_bin_form(a, b1, c1)
    return (a, b, c)

if __name__ == "__main__":
    print("example: check if forms (2,0,3) and (2,4,5) are equivalent")
    print(reduce_bin_form(2, 0, 3) == reduce_bin_form(2, 4, 5))