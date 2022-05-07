#!/usr/bin/env python

from dataclasses import dataclass

@dataclass(frozen=True, eq=True)
class Point:
    """
    Simple Point class with two coordinates x,y

    Attributes:
        x: int
        y: int
    """

    x: int = 0
    y: int = 0


class EllipticCurve:
    """
    Basic (Naive) implementation of the Elliptic Curves over the Field of Integers mod p (p for prime)

    Attributes
    ----------

        O: Point(None, None) 
            The point at infinity
        a, b, p: int
            Parameters of the EC: y^2 = x^3 + a*x + b mod p

    Methods: (incomplete, only the important ones)
    ----------

        __init__(a: int, b: int, p:int) -> None
            a, b: The coefficients for the poly
               p: Z/pZ
        add(P: Point, Q: Point) -> Point
            If two valid Point's are provided, returns the result of the Addition of `P`, `Q` on the defined EC 
        double(P: Point) -> Point
            If a valid `P` is a valid point, returns 2*P in the defined EC
        bruteforce_points() -> list
            Brutforce search all the points on the EC
        mod_add(P: Point, scalar: int) -> Point
            Double&Add algorithm
    """

    # Im. point at infinity
    O = Point(None, None)

    def __init__(self, a: int, b: int, p: int) -> None:
        self.__p = p
        self.__a = a
        self.__b = b

    def __repr__(self) -> str:
        return f'EllipticCurve(a={self.__a}, b={self.__b}, p={self.__p})'

    def __str__(self) -> str:
        """TeX repr. of the EC"""
        return f"y^2 \\heq x^3 + {self.__a: 3}\\cdot x + {self.__b: 3} \\mod {self.__p}"

    def double(self, P: Point):
        """Return the point `P` 'added' to itself"""

        if not self.valid(P):
            raise Exception(f"Given Point {P} is not on the EC {repr(self)}")

        if P == EllipticCurve.O:
            return EllipticCurve.O

        slope = ((3 * P.x**2 + self.__a) * self.inv_mod(2 * P.y)) % self.__p

        x = (slope**2 - 2 * P.x) % self.__p
        y = (slope * (P.x - x) - P.y) % self.__p

        return Point(x, y)

    def add(self, P: Point, Q: Point = None):
        """Return the 'sum' of two Points `P`, `Q`"""

        if not self.valid(P):
            raise Exception(f"Given Point {P} is not on the EC {repr(self)}")
        if not self.valid(Q):
            raise Exception(f"Given Point {Q} is not on the EC {repr(self)}")

        if P == Q:
            return self.double(P)
        if P == EllipticCurve.O:
            return Q
        elif Q == EllipticCurve.O:
            return P
        elif Q == self.inv(P):
            return EllipticCurve.O

        slope = ((Q.y - P.y) * self.inv_mod(Q.x - P.x)) % self.__p

        x = (slope**2 - (P.x + Q.x)) % self.__p
        y = (slope * (P.x - x) - P.y) % self.__p

        return Point(x, y)

    # ----------------------------------- Utils ---------------------------------- #
    def bruteforce_points(self) -> list:
        points = []
        for i in range(self.__p):
            for j in range(self.__p):
                P = Point(i, j)
                if self.valid(P):
                    points.append(P)
        return points

    def inv(self, P) -> Point:
        """Get the inverse of the point P on the defined EC"""

        if P == EllipticCurve.O:
            return P

        return Point(P.x, (-P.y) % self.__p)

    def valid(self, P: Point) -> bool:
        """Determine whether the the given point `P` lies on the curve"""

        if EllipticCurve.O == P:
            return True

        left = P.y**2 % self.__p
        right = (P.x**3 + (self.__a * P.x) + self.__b) % self.__p

        return left == right

    def crypto_check(self) -> bool:
        """Proofs, whether the parameters satisfy the property 4*a^3 + 27*b^2 != 0 mod p"""
        to_check = (4 * self.__a ** 3 + 27 * self.__b ** 2) % self.__p

        return (to_check, to_check != 0)
    
    def mod_add(self, P: Point, scalar: int) -> Point:
        """Square&Multiply, only it's Double&Add"""
        
        if not self.valid(P):
            raise Exception(f"Point {P} is not on the EC {repr(self)}")

        if scalar & 1:
            T = P
        else:
            T = EllipticCurve.O

        while scalar > 1:
            scalar >>= 1
            T = self.double(T)
            print(f'Double: {scalar} -> {T}')
            
            if (scalar & 1):
                T = self.add(T, P)
                print(f'Add: {scalar} -> {T}')

        return T
    
    def mod_pow(self, base, exp, mod=None) -> int:
        """Modular Exponentiation"""

        if mod is None:
            mod = self.__p

        if exp & 1:
            result = base
        else:
            result = 1

        while exp:
            exp >>= 1
            base = (base * base) % mod

            if (exp & 1):
                result = (result * base) % mod

        return result

    def inv_mod_pow(self, x: int) -> int:
        """Compute an inverse for x modulo `self.__p` (Little Fermat's Theo)"""

        if x % self.__p == 0:
            raise ZeroDivisionError(
                f"No Inverse for {x} in the Field of Integers mod {self.__p}")

        return self.mod_pow(x, self.__p-2, self.__p)

    def inv_mod(self, x: int) -> int:
        """Compute an inverse for x modulo `self.__p` (EEA)"""

        if x % self.__p == 0:
            raise ZeroDivisionError(
                f"No Inverse for {x} in the Field of Integers mod {self.__p}")

        return self.eea(self.__p, x)

    def eea(self, a: int, b: int):
        """Extended Euclidean Algorithm, a >= b"""
        if b < 0:
            b += a
        temp = a
        inverse = 0
        y = 1
        while b:
            q = a // b
            y, inverse = inverse - q * y, y
            a, b = b, a - q * b
        if inverse < 0:
            return inverse + temp
        return inverse

if __name__ == "__main__":
    
    pass
    
    """ Aufgabe 3, Teil (b) """

    """
        Point(x=278, y=796)  
        Point(x=805, y=129)  
        Point(x=None, y=None)
        Point(x=1088, y=390) 
    """
    
    """
    E = EllipticCurve(41, 297, 1217)

    args = \
        [
            (Point(*t[0]), Point(*t[1])) for t in
            [
                ((4, 1094), (4, 1094)),
                ((504, 1134), (526, 322)),
                ((271, 975), (271, 242)),
                ((1088, 390), (None, None))
            ]
        ]

    for arg in args:
        print(E.add(*arg))
    """

    # ------------ For test purposes: Results in `elliptic_curve.log` ------------ #
    # from logging import basicConfig as logconfig, INFO, info as loginfo
    # from itertools import product

    # def sep(len=31):
    #     return ' '.join('-' for _ in range(len))

    # assignment_format = f'{sep()}%(message)s'
    # logconfig(filename=f'.{__file__.split(".")[1]}.log', filemode='+w', level=INFO,
    #           format=assignment_format, datefmt='%b/%d/%Y %H:%M:%S')

    # sols_ec_17 = \
    #     [
    #         (None, None),
    #         (0, 6),
    #         (0, 11),
    #         (3, 1),
    #         (3, 16),
    #         (5, 1),
    #         (5, 16),
    #         (6, 3),
    #         (6, 14),
    #         (7, 6),
    #         (7, 11),
    #         (9, 1),
    #         (9, 16),
    #         (10, 6),
    #         (10, 11),
    #         (13, 7),
    #         (13, 10),
    #         (16, 4),
    #         (16, 13)
    #     ]

    # points = [Point(*sol) for sol in sols_ec_17]
    # E = EllipticCurve(*(2, 2), 17)

    # loginfo(f'\nAddition of Points on the Elliptic Curve\n{repr(E)} over the Field of Integers mod 17\nPoint(None, None) is the Im. point at infinity')
    # for (P, Q) in product(points, points):
    #     loginfo(f'\n{P} + {Q} = {E.add(P, Q)}')
    # ---------------------------------------------------------------------------- #
