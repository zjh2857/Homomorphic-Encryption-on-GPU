import tenseal as ts
import numpy as np

def gencontext():
    context = ts.context(ts.SCHEME_TYPE.CKKS, 8192, coeff_mod_bit_sizes=[22 ,21, 21, 21, 21, 21, 21, 21, 21, 21])
    context.global_scale = pow(2, 21)
    context.generate_galois_keys()
    return context

def encrypt(context, np_tensor):
    return ts.ckks_tensor(context, np_tensor)

def decrypt(enc_tensor):
    return np.array(enc_tensor.decrypt().tolist())

def bootstrap(context, tensor):
    # To refresh a tensor with exhausted depth. 
    # Here, bootstrap = enc(dec())
    tmp = decrypt(tensor)
    return encrypt(context, tmp)

if __name__ == "__main__":

    a = np.array([[1.,2.,3.,4.], [1.,2.,5.,4.]])
    context = gencontext()
    enc_a = encrypt(context, a)
    enc_at = encrypt(context, a.T)
    enc_b = encrypt(context, a)
    res = enc_a + enc_b
    res = enc_a - enc_b
    res = enc_a * enc_b
    res = enc_a @ enc_at
    print(decrypt(res))
