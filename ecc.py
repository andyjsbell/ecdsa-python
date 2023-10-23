from random import SystemRandom # cryptographic random byte generator

rand=SystemRandom() # create strong random number generator

# Convert a string with hex digits, colons, and whitespace to a long integer
def hex2int(hexString):
    return int("".join(hexString.replace(":","").split()),16)

# Half the extended Euclidean algorithm:
#    Computes   gcd(a,b) = a*x + b*y
#    Returns only gcd, x (not y)
# From http://rosettacode.org/wiki/Modular_inverse#Python
def half_extended_gcd(aa, bb):
    lastrem, rem = abs(aa), abs(bb)
    x, lastx = 0, 1
    while rem:
        lastrem, (quotient, rem) = rem, divmod(lastrem, rem)
        x, lastx = lastx - quotient*x, x
    return lastrem, lastx

# Modular inverse: compute the multiplicative inverse i of a mod m:
#     i*a = a*i = 1 mod m
def modular_inverse(a, m):
    g, x = half_extended_gcd(a, m)
    if g != 1:
        raise ValueError
    return x % m


# An elliptic curve has these fields:
#   p: the prime used to mod all coordinates
#   a: linear part of curve: y^2 = x^3 + ax + b
#   b: constant part of curve
#   G: a curve point (G.x,G.y) used as a "generator"
#   n: the order of the generator
class ECcurve:
    def __init__(self):
        return

    # Prime field multiplication: return a*b mod p
    def field_mul(self,a,b):
        return (a*b)%self.p

    # Prime field division: return num/den mod p
    def field_div(self,num,den):
        inverse_den=modular_inverse(den%self.p,self.p)
        return self.field_mul(num%self.p,inverse_den)

    # Prime field exponentiation: raise num to power mod p
    def field_exp(self,num,power):
        return pow(num%self.p,power,self.p)

    # Return the special identity point
    #   We pick x=p, y=0
    def identity(self):
        return ECpoint(self,self.p,0)

    # Return true if point Q lies on our curve
    def touches(self,Q):
        y2=self.field_exp(Q.y,2)
        x3ab=(self.field_mul((Q.x*Q.x)%self.p+self.a,Q.x)+self.b)%self.p
        return y2==x3ab

    # Return the slope of the tangent of this curve at point Q
    def tangent(self,Q):
        return self.field_div(Q.x*Q.x*3+self.a,Q.y*2)

    # Return the (x,y) point where this line intersects our curve
    #  Q1 and Q2 are two points on the line of slope m
    def line_intersect(self,Q1,Q2,m):
        v=(Q1.y + self.p - (m*Q1.x)%self.p)%self.p
        x=(m*m + self.p-Q1.x + self.p-Q2.x)%self.p
        y=(self.p-(m*x)%self.p + self.p-v)%self.p
        return ECpoint(self,x,y)

    # Return a doubled version of this elliptic curve point
    def double(self,Q):
        if (Q.x==self.p): # doubling the identity
            return Q
        return self.line_intersect(Q,Q,self.tangent(Q))

    # Return the "sum" of these elliptic curve points
    def add(self,Q1,Q2):
        # Identity special cases
        if (Q1.x==self.p): # Q1 is identity
            return Q2
        if (Q2.x==self.p): # Q2 is identity
            return Q1

        # Equality special cases
        if (Q1.x==Q2.x):
            if (Q1.y==Q2.y): # adding point to itself
                return self.double(Q1)
            else: # vertical pair--result is the identity
                return self.identity()

        # Ordinary case
        m=self.field_div(Q1.y+self.p-Q2.y,Q1.x+self.p-Q2.x)
        return self.line_intersect(Q1,Q2,m)

    # "Multiply" this elliptic curve point Q by the integer m
    #    Often the point Q will be the generator G
    def mul(self,Q,m):
        R=self.identity() # return point
        while m!=0:  # binary multiply loop
            if m&1: # bit is set
                # print("  mul: adding Q to R =",R);
                R=self.add(R,Q)
            m=m>>1
            if (m!=0):
                # print("  mul: doubling Q =",Q);
                Q=self.double(Q)

        return R

# A point on an elliptic curve: (x,y)
class ECpoint:
    """A point on an elliptic curve (x,y)"""
    def __init__(self,curve, x,y):
        self.curve=curve
        self.x=x
        self.y=y
        if not x==curve.p and not curve.touches(self):
            print(" ECpoint left curve: ",x,",",y)

    # "Add" this point to another point on the same curve
    def add(self,Q2):
        return self.curve.add(self,Q2)

    # "Multiply" this point by a scalar
    def mul(self,m):
        return self.curve.mul(self,m)

    # Print this ECpoint
    def __str__(self):
        if (self.x==self.curve.p):
            return "identity_point"
        else:
            return "("+str(self.x)+", "+str(self.y)+")"


