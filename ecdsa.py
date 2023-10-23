import random
from ecc import ECcurve, hex2int, ECpoint
import hashlib


def generate_curve():
    # defaults to a secp256k1 curve
    secp256k1 = ECcurve()
    secp256k1.p = hex2int("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F")
    secp256k1.a = 0
    secp256k1.b = 7
    secp256k1.n = 115792089237316195423570985008687907852837564279074904382605163141518161494337

    secp256k1.G = ECpoint(curve=secp256k1,
                          x=55066263022277343669578718895168534326250603453777594175500187360389116729240,
                          y=32670510020758816978083085130507043184471273380659243275938904335757337482424
                          )
    return secp256k1


def generate_keypair():
    curve = generate_curve()

    # Generate key pair
    private_key = random.randint(0, curve.n)
    public_key = curve.mul(curve.G, private_key)

    return curve, public_key, private_key


def generate_hash(plain_message):
    return int.from_bytes(hashlib.sha256(plain_message).digest(), byteorder='big')


def sign_message(curve, plain_message, private_key):
    # Hash plain message
    h = generate_hash(plain_message)
    # Sign message with private key
    # Calculate the random point R = k * G and take its x-coordinate: r = R.x
    k = random.randint(0, curve.n)
    R = curve.mul(curve.G, k)
    r = R.x
    # Calculate the signature proof: s = k^−1 * (h+r*privKey)(mod n)
    s = pow(k, -1, curve.n) * (h + r * private_key) % curve.n

    # Return the signature
    return r, s


def verify_signature(curve, public_key, r, s, plain_message):
    # Hash the message
    h = generate_hash(plain_message)
    # Calculate the modular inverse of the signature proof s1 = s^-1 mod n
    s1 = pow(s, -1, curve.n)
    # Recover the random point used during the signing: R' = (h * s1) * G + (r * s1) * pubKey
    R = curve.add(curve.mul(curve.G, h * s1), curve.mul(public_key, r * s1))

    return R.x == r


if __name__ == '__main__':
    # Generate key pair
    (curve, public_key, private_key) = generate_keypair()

    # Generate a hash of our message
    plain_message = "all work and no play makes jack a dull boy".encode()
    r, s = sign_message(curve, plain_message, private_key)
    print(f"Message {plain_message} signed producing signature:\nr = {r}, s = {s}")
    print("Verifying signature....")
    assert verify_signature(curve, public_key, r, s, plain_message)
    print("Verified message with signature!")