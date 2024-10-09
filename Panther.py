# Created by Alexandre Mihet, Billal Medour, Jean-Paul Morcos Doueihy, Malak Lalami


# **************************************
# **************************************
# General functions
# **************************************
# **************************************

def hexa_to_binary(hexa_string):
    """
    Takes a string of hexadecimal numbers and returns its binary representation.

    Args:
        hexa_string: A string of hexadecimal numbers.
    
    Results:
        result (str): A string of binary numbers.
    """
    result = ""
    for char in hexa_string:
        result += bin(int(char, 16))[2:].zfill(4)

    return result

def binary_to_hexa(binary_num):
    """
    Takes a string or a list of binary numbers and returns its hexadecimal representation.

    Args:
        binary_num: A string or a list of binary numbers.
    
    Results:
        result (str): A string of hexadecimal numbers.
    """
    binary_num = ''.join(str(i) for i in binary_num)
    result = hex(int(binary_num, 2))[2:]

    return result

# **************************************
# **************************************
# Cipher implementation
# **************************************
# **************************************

def get_sr(state):
    """
    Takes the state as input and returns the outer part 'sr' as a list.

    Args:
        state: A string representing the state of 328 bits.
    
    Returns:
        list: A list of 64 characters representing the outer part 'sr' of the state.
    """
    sr = ""
    sr = state[60:76] + state[140:156] + state[224:240] + state[312:328]
    sr = list(sr)
    return sr

def set_sr(state, sr):
    """
    Takes the state and it's outer part 'sr' as inputs and returns the state after updating it.

    Args:
        state: A string representing the state of 328 bits.
        sr: A list of 64 character representing the outer part 'sr' of the state.
    
    Returns:
        state (str): the state of 328 bits after updating it's outer part 'sr'.
    """
    sr = ''.join(str(bit) for bit in sr)
    state = state[0:60] + sr[0:16] + state[76:140] + sr[16:32] + state[156:224] + sr[32:48] + state[240:312] + sr[48:] 
    return state

def toeplitz_multiplication(nibble_1, nibble_2, nibble_3, nibble_4):
    """
    Takes four nibbles (strings of four bits) as input and returns four other nibbles after performing a matrix multiplication.

    Args:
        nibble_1, nibble_2, nibble_3, nibble_4: each one of the inputs is a string of four bits representing a hexadecimal number.
    
    Returns:
        output_1, output_2, output_3, output_4 (str): each one of the outputs is a string of four bits representing a hexadecimal number after the matrix multiplication.

    P.S:\n
    - In the equations of the function, the inputs are referred as "l1, l2, l3, l4" and the outputs as "d1, d2, d3, d4".\n
    - Each nibble is considered a third degree polynomial, nibble = "e3e2e1e0" = e3 * X^3 + e2 * X^2 + e1 * X + e0\n
    - Calculations have been made and these are the equations needed:\n
        d1 = l1 + l2 + l3 * X + l4 * (X^3 + X^2)\n
        d2 = l1(X^2 + X) + l2 + l3 + l4 * X\n
        d3 = l1 + l2(X^2 + X) + l3 + l4\n
        d4 = l1(X^3 + X^2) + l2 + l3(X^2 + X) + l4\n
    having the following:\n
        X * l = (e3+e2)X^3 + e1 * X^2 + e0 * X + e3\n
        (X^3 + X^2) * l = e0 * X^3 + (e3+e0)X^2 + e2 * X + e1\n
        (X^2 + X) * l = e1 * X^3 + (e1+e0)X^2 + (e0+e3)X + e2\n
    """
    l1, l2, l3, l4 = nibble_1, nibble_2, nibble_3, nibble_4
    
    # d1 = l1 + l2 + l3 * X + l4 * (X^3 + X^2)
    d1_l1 = l1
    d1_l2 = l2
    d1_l3 = '{}{}{}{}'.format((int(l3[0], 2) + int(l3[1], 2)) % 2, l3[2], l3[3], l3[0])
    d1_l4 = '{}{}{}{}'.format(l4[3], (int(l4[0], 2) + int(l4[3], 2)) % 2, l4[1], l4[2])
    d1 = bin(int(d1_l1, 2) ^ int(d1_l2, 2) ^ int(d1_l3, 2) ^ int(d1_l4, 2))[2:]
    d1 = d1.zfill(4)

    # d2 = l1(X^2 + X) + l2 + l3 + l4 * X
    d2_l1 = '{}{}{}{}'.format(l1[2], (int(l1[2], 2) + int(l1[3], 2)) % 2, (int(l1[3], 2) + int(l1[0], 2)) % 2, l1[1])
    d2_l2 = l2
    d2_l3 = l3
    d2_l4 = '{}{}{}{}'.format((int(l4[0], 2) + int(l4[1], 2)) % 2, l4[2], l4[3], l4[0])
    d2 = bin(int(d2_l1, 2) ^ int(d2_l2, 2) ^ int(d2_l3, 2) ^ int(d2_l4, 2))[2:]
    d2 = d2.zfill(4)

    # d3 = l1 + l2(X^2 + X) + l3 + l4
    d3_l1 = l1
    d3_l2 = '{}{}{}{}'.format(l2[2], (int(l2[2], 2) + int(l2[3], 2)) % 2, (int(l2[3], 2) + int(l2[0], 2)) % 2, l2[1])
    d3_l3 = l3
    d3_l4 = l4
    d3 = bin(int(d3_l1, 2) ^ int(d3_l2, 2) ^ int(d3_l3, 2) ^ int(d3_l4, 2))[2:]
    d3 = d3.zfill(4)

    # d4 = l1(X^3 + X^2) + l2 + l3(X^2 + X) + l4
    d4_l1 = '{}{}{}{}'.format(l1[3], (int(l1[0], 2) + int(l1[3], 2)) % 2, l1[1], l1[2])
    d4_l2 = l2
    d4_l3 = '{}{}{}{}'.format(l3[2], (int(l3[2], 2) + int(l3[3], 2)) % 2, (int(l3[3], 2) + int(l3[0], 2)) % 2, l3[1])
    d4_l4 = l4
    d4 = bin(int(d4_l1, 2) ^ int(d4_l2, 2) ^ int(d4_l3, 2) ^ int(d4_l4, 2))[2:]
    d4 = d4.zfill(4)

    return d1, d2, d3, d4

def field_multiplication(a, b):
    """
    Takes two nibbles (strings of four bits) as input and returns a nibble representing the multiplication of the inputs.

    Args:
        a, b: each one of the inputs is a string of four bits representing a hexadecimal number.
    
    Returns:
        result (str): a string of four bits representing the multiplication of the inputs.

    P.S:\n
    - Each nibble is considered a third degree polynomial. These are the equations used to form "dict_a" and dict_b":\n
        a = "a3a2a1a0" = a3 * X^3 + a2 * X^2 + a1 * X + a0\n
        b = "b3b2b1b0" = b3 * X^3 + b2 * X^2 + b1 * X + b0\n

    - The following equation is used to form "result_dict":\n
    a * b = a3b3X6 + a3b2X5 + a3b1X4 + a3b0X3 + a2b3X5 + a2b2X4 + a2b1X3 + a2b0X2 + a1b3X4 + a1b2X3 + a1b1X2 + a1b0X + a0b3X3 + a0b2X2 + a0b1X + a0b0\n
        = X6(a3b3) + X5(a3b2 + a2b3) + X4(a3b1 + a2b2 + a1b3) + X3(a3b0 + a2b1 + a1b2 + a0b3) + X2(a2b0 + a1b1 + a0b2) + X1(a1b0 + a0b1) + X0(a0b0)\n
    having the following (used in the if clauses):\n
        X4 = X3 + 1\n
        X5 = X3 + X + 1\n
        X6 = X3 + X2 + X + 1\n
    """
    a3, a2, a1, a0 = a[0], a[1], a[2], a[3]
    b3, b2, b1, b0 = b[0], b[1], b[2], b[3]

    dict_a = {0 : int(a0),
              1 : int(a1),
              2 : int(a2),
              3 : int(a3)}
    dict_b = {0 : int(b0),
              1 : int(b1),
              2 : int(b2),
              3 : int(b3)}
    result_dict = {0 : (dict_a[0] * dict_b[0]),
                   1 : ((dict_a[1] * dict_b[0] + dict_a[0] * dict_b[1])) % 2,
                   2 : ((dict_a[2] * dict_b[0] + dict_a[1] * dict_b[1] + dict_a[0] * dict_b[2]) % 2),
                   3 : ((dict_a[3] * dict_b[0] + dict_a[2] * dict_b[1] + dict_a[1] * dict_b[2] + dict_a[0] * dict_b[3]) % 2),
                   4 : ((dict_a[3] * dict_b[1] + dict_a[2] * dict_b[2] + dict_a[1] * dict_b[3])) % 2,
                   5 : ((dict_a[3] * dict_b[2] + dict_a[2] * dict_b[3])) % 2,
                   6 : (dict_a[3] * dict_b[3])}
    if result_dict[4] == 1:
        result_dict[4] = 0
        result_dict[3] = (result_dict[3] + 1) % 2
        result_dict[0] = (result_dict[0] + 1) % 2
    if result_dict[5] == 1:
        result_dict[5] = 0
        result_dict[3] = (result_dict[3] + 1) % 2
        result_dict[1] = (result_dict[1] + 1) % 2
        result_dict[0] = (result_dict[0] + 1) % 2
    if result_dict[6] == 1:
        result_dict[6] = 0
        result_dict[3] = (result_dict[3] + 1) % 2
        result_dict[2] = (result_dict[2] + 1) % 2
        result_dict[1] = (result_dict[1] + 1) % 2
        result_dict[0] = (result_dict[0] + 1) % 2

    result = ''.join(str(result_dict[i]) for i in range(3, -1, -1))

    return result

def F(state, n):
    """
    Takes the state and the number of iterations this function does 'n' and returns the updated state.

    Args:
        state: A string representing the state of 328 bits.
        n: An integer representing the number of iterations this function does.
    
    Returns:
        state (str): the state of 328 bits after updating it.
    """
    state_len = 328
    for _ in range(n):
        state_nibbles = [state[i:i + 4] for i in range(0, state_len, 4)]

        P = state_nibbles[0:19]
        Q = state_nibbles[19:39]
        R = state_nibbles[39:60]
        S = state_nibbles[60:82]
        
        # Calculate feedback polynomial
        fp = int(P[0], 2) ^ int(P[7], 2) ^ int(P[10], 2) ^ (int(field_multiplication(P[6], P[18]), 2))
        fq = int(Q[0], 2) ^ int(Q[4], 2) ^ int(Q[6], 2) ^ int(Q[7], 2) ^ int(Q[15], 2) ^ (int(field_multiplication(Q[3], Q[7]), 2))
        fr = int(R[0], 2) ^ int(R[1], 2) ^ int(R[15], 2)  ^ int(R[17], 2)  ^ int(R[19], 2) ^ (int(field_multiplication(R[13], R[15]), 2))
        fs = int(S[0], 2) ^ int(S[1], 2) ^ (int(field_multiplication(S[4], S[10]), 2)) ^ (int(field_multiplication(S[11], S[18]), 2))

        # Calculate interconncection polynomial
        gp = int(Q[9], 2) ^ int(R[10], 2) ^ int(S[12], 2)
        gq = int(P[4], 2) ^ int(R[2], 2) ^ int(S[5], 2)
        gr = int(P[12], 2)  ^ int(Q[11], 2)  ^ int(S[16], 2)
        gs = int(P[16], 2)  ^ int(Q[17], 2)  ^ int(R[2], 2) 
        
        # Round constant selection
        rc1 = 7
        rc2 = 9
        rc3 = 11
        rc4 = 13
        # Create inputs for filter function
        l1 = bin(fp ^ gp ^ rc1)[2:].zfill(4)
        l2 = bin(fq ^ gq ^ rc2)[2:].zfill(4)
        l3 = bin(fr ^ gr ^ rc3)[2:].zfill(4)
        l4 = bin(fs ^ gs ^ rc4)[2:].zfill(4)
        # Toeplitz multiplication
        d1, d2, d3, d4 = toeplitz_multiplication(l1, l2, l3, l4)
        # S-box followed by Toeplitz multiplication
        t1, t2, t3, t4 = toeplitz_multiplication(Substitution(d1), Substitution(d2), Substitution(d3), Substitution(d4))
        # Shift the Registers and add nibbles
        state = state[4:76] + t1 + state[80:156] + t2 + state[160:240] + t3 + state[244:328] + t4        
    
    return state

def Substitution(nibble):
    """
    Takes a string of four bits (nibble) as input and returns the corresponding substitution.

    Args:
        nibble: A string of four bits.

    Returns:
        result (str): A string of four bits representing the corresponding substitution of the input nibble.
    """
    S = {
        '0000':'1100',
        '0001':'0101',
        '0010':'0110',
        '0011':'1011',
        '0100':'1001',
        '0101':'0000',
        '0110':'1010',
        '0111':'1101',
        '1000':'0011',
        '1001':'1110',
        '1010':'1111',
        '1011':'1000',
        '1100':'0100',
        '1101':'0111',
        '1110':'0001',
        '1111':'0010'}

    result = S[nibble]
    return result

def Initialization(key, IV):
    """
    Takes the key and the initial value 'IV' as inputs and returns the state.

    Args:
        key: A string representing the key of 128 bits.
        IV: A string representing the initial value of 128 bits.

    Returns:
        state (str): the initial state containing "the key || the IV || the first 64 bits of the complement of the key || seven 1's || a 0".
    """
    state = ""
    for i in range(0, 128):
        state += key[i]

    for i in range(0, 128):
        state += IV[i]

    for i in range(0, 64):
        state += str(1 - int(key[i]))

    for i in range(0, 7):
        state += "1"

    state += "0"
    state = F(state, 92)

    return state

def ADProcessing(state, AD):
    """
    Takes the state and the associated data 'AD' as inputs and returns the state.

    Args:
        state: A string representing the state of 328 bits.
        AD: A string representing the associated data. Its length can vary depending on the input given but will be padded to be a multiple of 64.

    Returns:
        state (str): the state of 328 bits after the absorbtion of the associated data.
    """
    # Padding if the length of the associated data is not a multiple of 64
    AD_len = len(AD)
    if AD_len % 64 != 0:
        pad = 64 - (AD_len % 64)-1
        AD = AD + "1" + "0" * pad
    
    # Divide the string 'AD' into blocks of 64 characters and gather the blocks in a list 'AD_blocks'
    # AD1,AD2,....,ADk = AD
    block_len = 64
    AD_blocks = []
    nb_blocks = 0
    for i in range(0, AD_len, block_len):
        AD_blocks.append(AD[i:i + block_len])
        nb_blocks += 1

    # Absorb the associated data block by block into the state
    for i in range(nb_blocks):
        sr = get_sr(state)
        for j in range(64):
            sr[j] = int(sr[j], 2) ^ int(AD_blocks[i][j], 2)
        state = set_sr(state, sr)
        state = F(state, 4)
    
    return state

def PlaintextProcessing(state, PT):
    """
    Takes the state and the plaintext 'PT' as inputs and returns the state and the ciphertext.

    Args:
        state: A string representing the state of 328 bits.
        PT: A string representing the plaintext. Its length can vary depending on the input given but will be padded to be a multiple of 64.

    Returns:
        state (str): the state of 328 bits after the encryption of the plaintext.
        CT (list): the ciphertext blocks gathered in a list.
    """
    # Padding if the length of the plaintext is not a multiple of 64
    PT_len = len(PT)
    if PT_len % 64 != 0:
        pad = 64 -(PT_len % 64) - 1
        PT = PT + "1" + "0" * pad

    # Divide the string 'PT' into blocks of 64 characters and gather the blocks in a list 'PT_blocks'
    # PT1,PT2, ··· ,PTn = PT
    block_len = 64
    PT_blocks = []
    nb_blocks = 0
    for i in range(0, PT_len, block_len):
        PT_blocks.append(PT[i:i + block_len])
        nb_blocks += 1

    # Get the corresponding ciphertext for each block of plaintext and gather them in a list 'CT'
    CT = []
    sr = get_sr(state)
    for i in range(nb_blocks):
        for j in range(64):
            sr[j] = int(sr[j], 2) ^ int(PT_blocks[i][j], 2)
            CT.append(sr[j])
        state = set_sr(state, sr)
        if i+1 < nb_blocks:
            state = F(state, 4)
            sr = get_sr(state)
        
    return state, CT

def CiphertextProcessing(state, CT):
    """
    Takes the state and the ciphertext 'CT' as inputs and returns the state and the plaintext.

    Args:
        state: A string representing the state of 328 bits.
        CT: A string representing the ciphertext. Its length can vary depending on the plaintext but is always a multiple of 64.

    Returns:
        state (str): the state of 328 bits after the decryption of the ciphertext.
        PT (list): the plaintext blocks gathered in a list.
    """
    # Divide the string 'CT' into blocks of 64 characters and gather the blocks in a list 'CT_blocks'
    block_len = 64
    CT_len = len(CT)
    CT_blocks = []
    nb_blocks = 0 
    for i in range(0, CT_len, block_len):
        CT_blocks.append(CT[i:i + block_len])
        nb_blocks += 1
    
    # Get the corresponding plaintext for each block of ciphertext and gather them in a list 'PT'
    PT = []
    sr = get_sr(state)
    for i in range(nb_blocks):
        for j in range(64):
            PT.append(int(sr[j]) ^ int(CT_blocks[i][j]))
            sr[j] = CT_blocks[i][j]
        state = set_sr(state, sr)
        if i+1 < nb_blocks:
            state = F(state, 4)
            sr = get_sr(state)
    
    return state, PT

def Finalization(state, hashlen):
    """
    Takes the state and a hashlen as inputs and returns the tag.

    Args:
        state: A string representing the state of 328 bits.
        hashlen: An integer representing the length of the tag in a binary representation.

    Returns:
        tag (list): the tag blocks gathered in a list.
    """
    state = F(state, 92)

    # Treating two cases: hashlen is a multiple of 64 (if clause), or not (else clause)
    tag = []
    sr = get_sr(state)
    if hashlen % 64 == 0:
        for _ in range(0, hashlen // 64):
            for j in range(0, 64):
                tag.append(sr[j])
            state = F(state, 4)
            sr = get_sr(state)
    else:
        for _ in range(0, hashlen // 64):
            for j in range(0, 64):
                tag.append(sr[j])
            state = F(state, 4)
            sr = get_sr(state)

        for j in range(0, hashlen % 64):
            tag.append(sr[j])

    return tag

def Encryption(key, IV, AD, PT, hashlen):
    """
    The encryption algorithm.

    Args:
        key: A string representing the key of 128 bits.
        IV: A string representing the initial value of 128 bits.
        AD: A string representing the associated data.
        PT: A string representing the plaintext.
        hashlen: An integer representing the length of the tag in a binary representation.

    Returns:
        CT (list): the ciphertext blocks gathered in a list.
        tag (list): the tag blocks gathered in a list.
    """
    state = 0
    state = Initialization(key, IV)
    state = ADProcessing(state, AD)
    state, CT = PlaintextProcessing(state, PT)
    tag = Finalization(state, hashlen)
    return CT, tag

def Decryption(key, IV, AD, CT, tag, hashlen):
    """
    The decryption algorithm.

    Args:
        key: A string representing the key of 128 bits.
        IV: A string representing the initial value of 128 bits.
        AD: A string representing the associated data.
        CT: A string representing the ciphertext.
        tag: A string representing the tag obtained after encrypting.
        hashlen: An integer representing the length of the tag in a binary representation.

    Returns:
        PT (list): the plaintext blocks gathered in a list.
        "Error": if the encryption tag does not match the decryption tag.
    """
    state = Initialization(key, IV)
    state = ADProcessing(state, AD)
    state, PT = CiphertextProcessing(state, CT)
    decTag = binary_to_hexa(Finalization(state, hashlen))
    tag = binary_to_hexa(tag)
    if tag == decTag:
        return PT
    else:
        return "Error"

# **************************************
# **************************************
# Attacks on cipher
# **************************************
# **************************************

# Inverse functions
def Inv_substitution(nibble):
    """
    Takes a string of four bits (nibble) as input and returns the corresponding substitution.

    Args:
        nibble: A string of four bits.

    Returns:
        result (str): A string of four bits representing the corresponding substitution of the input nibble.
    """
    S = {
        '0000': '0101',
        '0001':'1110',
        '0010':'1111',
        '0011':'1000',
        '0100':'1100',
        '0101':'0001',
        '0110':'0010',
        '0111':'1101',
        '1000':'1011',
        '1001':'0100',
        '1010':'0110',
        '1011':'0011',
        '1100':'0000',
        '1101':'0111',
        '1110':'1001',
        '1111':'1010'}

    result = S[nibble]
    return result

def Inv_toeplitz_multiplication(d1, d2, d3, d4):
    """
    Takes four nibbles (strings of four bits) as input and returns four other nibbles after performing a matrix multiplication.

    Args:
        nibble_1, nibble_2, nibble_3, nibble_4: each one of the inputs is a string of four bits representing a hexadecimal number.
    
    Returns:
        output_1, output_2, output_3, output_4 (str): each one of the outputs is a string of four bits representing a hexadecimal number after the matrix multiplication.

    P.S:\n
    - In the equations of the function, the inputs are referred as "d1, d2, d3, d4" and the outputs as "l1, l2, l3, l4".\n
    - Each nibble is considered a third degree polynomial, nibble = "e3e2e1e0" = e3 * X^3 + e2 * X^2 + e1 * X + e0\n
    - Calculations have been made and these are the equations needed:\n
        l1 = (X^3 + X^2 + X)d1 + (X^3 + X^2)d2 + (X^3 + X)d3 + (X^3 + X^2)d4\n
        l2 = (X^2 + X)d1 + (X^3 + X^2 + X)d2 + (X^3 + X^2)d3 + (X^3 + X)d4\n
        l3 = (X^2 + 1)d1 + (X^2 + X)d2 + (X^3 + X^2 + X)d3 + (X^3 + X^2)d4\n
        l4 = (X^2 + X)d1 + (X^2 + 1)d2 + (X^2 + X)d3 + (X^3 + X^2 + X)d4\n
    having the following:\n
        (X^2 + 1)d = (e2 + e1)X^3 + (e2 + e0)X^2 + (e3 + e1)X + (e3 + e2 + e0)\n
        (X^3 + X)d = (e1 + e0)X^3 + (e3 + e1)X^2 + (e3 + e2 + e0)X + (e2 + e1)\n
        (X^3 + X^2 + X)d = (e3 + e2 + e0)X^3 + (e3 + e1 + e0)X^2 + (e2 + e0)X + (e3 + e1)\n
        (X^3 + X^2)d = e0*X^3 + (e3 + e0)X^2 + e2*X + e1\n
        (X^2 + X)d = e1*X^3 + (e1 + e0)X^2 + (e3 + e0)X + e2\n
    """

    # l1 = (X^3 + X^2 + X)d1 + (X^3 + X^2)d2 + (X^3 + X)d3 + (X^3 + X^2)d4
    l1_d1 = '{}{}{}{}'.format((int(d1[0], 2) + int(d1[1], 2) + int(d1[3], 2)) % 2,
                           (int(d1[0], 2) + int(d1[2], 2) + int(d1[3], 2)) % 2,
                           (int(d1[1], 2) + int(d1[3], 2)) % 2,
                           (int(d1[0], 2) + int(d1[2], 2)) % 2)

    l1_d2 = '{}{}{}{}'.format(d2[3], (int(d2[0], 2) + int(d2[3], 2)) % 2, d2[1], d2[2])

    l1_d3 = '{}{}{}{}'.format((int(d3[2], 2) + int(d3[3], 2)) % 2,
                           (int(d3[0], 2) + int(d3[2], 2)) % 2,
                           (int(d3[0], 2) + int(d3[1], 2) + int(d3[3], 2)) % 2,
                           (int(d3[1], 2) + int(d3[2], 2)) % 2)
    
    l1_d4 = '{}{}{}{}'.format(d4[3], (int(d4[0], 2) + int(d4[3], 2)) % 2, d4[1], d4[2])
    
    l1 = bin(int(l1_d1, 2) ^ int(l1_d2, 2) ^ int(l1_d3, 2) ^ int(l1_d4, 2))[2:]
    l1 = l1.zfill(4)

    # l2 = (X^2 + X)d1 + (X^3 + X^2 + X)d2 + (X^3 + X^2)d3 + (X^3 + X)d4
    l2_d1 = '{}{}{}{}'.format(d1[2], (int(d1[2], 2) + int(d1[3], 2)) % 2, (int(d1[0], 2) + int(d1[3], 2)) % 2, d1[1])

    l2_d2 = '{}{}{}{}'.format((int(d2[0], 2) + int(d2[1], 2) + int(d2[3], 2)) % 2,
                           (int(d2[0], 2) + int(d2[2], 2) + int(d2[3], 2)) % 2,
                           (int(d2[1], 2) + int(d2[3], 2)) % 2,
                           (int(d2[0], 2) + int(d2[2], 2)) % 2)

    l2_d3 = '{}{}{}{}'.format(d3[3], (int(d3[0], 2) + int(d3[3], 2)) % 2, d3[1], d3[2])
    
    l2_d4 = '{}{}{}{}'.format((int(d4[2], 2) + int(d4[3], 2)) % 2,
                           (int(d4[0], 2) + int(d4[2], 2)) % 2,
                           (int(d4[0], 2) + int(d4[1], 2) + int(d4[3], 2)) % 2,
                           (int(d4[1], 2) + int(d4[2], 2)) % 2)
    
    l2 = bin(int(l2_d1, 2) ^ int(l2_d2, 2) ^ int(l2_d3, 2) ^ int(l2_d4, 2))[2:]
    l2 = l2.zfill(4)

    # l3 = (X^2 + 1)d1 + (X^2 + X)d2 + (X^3 + X^2 + X)d3 + (X^3 + X^2)d4
    l3_d1 = '{}{}{}{}'.format((int(d1[1], 2) + int(d1[2], 2)) % 2,
                           (int(d1[1], 2) + int(d1[3], 2)) % 2,
                           (int(d1[0], 2) + int(d1[2], 2)) % 2,
                           (int(d1[0], 2) + int(d1[1], 2) + int(d1[3], 2)) % 2)

    l3_d2 = '{}{}{}{}'.format(d2[2], (int(d2[2], 2) + int(d2[3], 2)) % 2, (int(d2[0], 2) + int(d2[3], 2)) % 2, d2[1])

    l3_d3 = '{}{}{}{}'.format((int(d3[0], 2) + int(d3[1], 2) + int(d3[3], 2)) % 2,
                           (int(d3[0], 2) + int(d3[2], 2) + int(d3[3], 2)) % 2,
                           (int(d3[1], 2) + int(d3[3], 2)) % 2,
                           (int(d3[0], 2) + int(d3[2], 2)) % 2)
    
    l3_d4 = '{}{}{}{}'.format(d4[3], (int(d4[0], 2) + int(d4[3], 2)) % 2, d4[1], d4[2])
    
    l3 = bin(int(l3_d1, 2) ^ int(l3_d2, 2) ^ int(l3_d3, 2) ^ int(l3_d4, 2))[2:]
    l3 = l3.zfill(4)
    
    # l4 = (X^2 + X)d1 + (X^2 + 1)d2 + (X^2 + X)d3 + (X^3 + X^2 + X)d4
    l4_d1 = '{}{}{}{}'.format(d1[2], (int(d1[2], 2) + int(d1[3], 2)) % 2, (int(d1[0], 2) + int(d1[3], 2)) % 2, d1[1])

    l4_d2 = '{}{}{}{}'.format((int(d2[1], 2) + int(d2[2], 2)) % 2,
                           (int(d2[1], 2) + int(d2[3], 2)) % 2,
                           (int(d2[0], 2) + int(d2[2], 2)) % 2,
                           (int(d2[0], 2) + int(d2[1], 2) + int(d2[3], 2)) % 2)

    l4_d3 = '{}{}{}{}'.format(d3[2], (int(d3[2], 2) + int(d3[3], 2)) % 2, (int(d3[0], 2) + int(d3[3], 2)) % 2, d3[1])
    
    l4_d4 = '{}{}{}{}'.format((int(d4[0], 2) + int(d4[1], 2) + int(d4[3], 2)) % 2,
                           (int(d4[0], 2) + int(d4[2], 2) + int(d4[3], 2)) % 2,
                           (int(d4[1], 2) + int(d4[3], 2)) % 2,
                           (int(d4[0], 2) + int(d4[2], 2)) % 2)
    
    l4 = bin(int(l4_d1, 2) ^ int(l4_d2, 2) ^ int(l4_d3, 2) ^ int(l4_d4, 2))[2:]
    l4 = l4.zfill(4)

    return l1, l2, l3, l4
    
def Inv_f(state, n):
    """
    Takes the state and the number of iterations this function does 'n' and returns the updated state.

    Args:
        state: A string representing the state of 328 bits.
        n: An integer representing the number of iterations this function does.
    
    Returns:
        state (str): the state of 328 bits after updating it.
    """
    state_len = 328
    for _ in range(n):
        state_nibbles = [state[i:i + 4] for i in range(0, state_len, 4)]

        P = state_nibbles[0:19]
        Q = state_nibbles[19:39]
        R = state_nibbles[39:60]
        S = state_nibbles[60:82]

        # Inverse Toeplitz multiplication and inverse substitution
        t1, t2, t3, t4 = P[18], Q[19], R[20], S[21]  
        e1, e2, e3, e4 = Inv_toeplitz_multiplication(t1, t2, t3, t4)
        d1, d2, d3, d4 = Inv_substitution(e1), Inv_substitution(e2), Inv_substitution(e3), Inv_substitution(e4)
        l1, l2, l3, l4 = Inv_toeplitz_multiplication(d1, d2, d3, d4)

        # Round constant selection
        rc1 = 7
        rc2 = 9
        rc3 = 11
        rc4 = 13

        # Calculate P0, Q0, R0 and S0
        P0 = int(l1, 2) ^ int(P[6], 2) ^ int(P[9], 2)^ (int(field_multiplication(P[5], P[17]), 2)) ^ int(Q[8], 2) ^ int(R[9], 2) ^ int(S[11], 2) ^ rc1
        Q0 = int(l2, 2) ^ int(Q[3], 2) ^ int(Q[5], 2) ^ int(Q[6], 2) ^ int(Q[14], 2) ^ (int(field_multiplication(Q[2], Q[6]), 2)) ^ int(P[3], 2) ^ int(R[1], 2) ^ int(S[4], 2) ^ rc2
        R0 = int(l3, 2) ^ int(R[0], 2) ^ int(R[14], 2) ^ int(R[16], 2) ^ int(R[18], 2) ^ (int(field_multiplication(R[12], R[14]), 2)) ^ int(P[11], 2) ^ int(Q[10], 2) ^ int(S[15], 2) ^ rc3
        S0 = int(l4, 2) ^ int(S[0], 2) ^ (int(field_multiplication(S[3], S[9]), 2)) ^ (int(field_multiplication(S[10], S[17]), 2)) ^ int(P[15], 2) ^ int(Q[16], 2) ^ int(R[1], 2) ^ rc4

        P0 = bin(P0)[2:].zfill(4)
        Q0 = bin(Q0)[2:].zfill(4)
        R0 = bin(R0)[2:].zfill(4)
        S0 = bin(S0)[2:].zfill(4)

        state = P0 + state[0:72] + Q0 + state[76:152] + R0 + state[156:236] + S0 + state[240:324]

    return state

# Key recovery attack
def Key_recovery_attack(PT, CT):
    """
    Takes the plaintext and the ciphertext as inputs and returns the key used for the encryption.

    Args:
        PT: A string representing the plaintext.
        CT: A string representing the ciphertext.

    Returns:
        recovered_key (str): A string representing the recovered key.
    
    P.S: In order for the attack to work, at least six 64-bit blocks of plaintext/ciphertext are needed.
    """
    # Recover the full state from the ciphertext
    tmp = ""
    for i in range(4):
        x = i*16
        tmp += CT[0+x:16+x] + CT[64+x:80+x] + CT [128+x:144+x] + CT[192+x:208+x] + CT[256+x:272+x] + CT[320+x:336+x]
    state = tmp[20:96] + tmp[112:192] + tmp[204:288] + tmp[296:384]

    # Divide the plaintext into blocks of 64 bits
    PT_blocks  = []
    for i in range(0, 384, 64):
        PT_blocks.append(PT[i:i+64])
    
    # XOR between the outer part of the state 'sr' and the plaintext blocks
    for i in range(5, 0, -1):
        sr = get_sr(state)
        for j in range (64):
            sr[j] = int(sr[j], 2) ^ int(PT_blocks[i][j], 2)
        state = set_sr(state, sr)
        state = Inv_f(state, 4)

    sr = get_sr(state)
    for j in range (64):
        sr[j] = int(sr[j], 2) ^ int(PT_blocks[0][j], 2)
    state = set_sr(state, sr)

    state = Inv_f(state, 92)

    recovered_key = state[0:128]

    return recovered_key

# Plaintext recovery attack
def Plaintext_recovery_attack(CT):
    """
    Takes the ciphertext as input and returns the recovered plaintext.

    Args:
        CT: A string representing the ciphertext.
    
    Returns:
        recovered_PT: A string representing the plaintext recovered.

    P.S: In this attack, only the blocks of plaintext after the sixth one are recovered.
    """
    # Recover the full state from the ciphertext
    tmp = ""
    for i in range(4):
        x = i*16
        tmp += CT[0+x:16+x] + CT[64+x:80+x] + CT [128+x:144+x] + CT[192+x:208+x] + CT[256+x:272+x] + CT[320+x:336+x]
    state = tmp[20:96] + tmp[112:192] + tmp[204:288] + tmp[296:384]

    # Divide the ciphertext into blocks of 64 bits
    CT_len = len(CT)
    nb_blocks = (CT_len//64)
    CT_blocks  = []
    for i in range(0, CT_len, 64):
        CT_blocks.append(CT[i:i+64])

    # Recover the plaintext blocks > 6; PT6...PTn
    recovered_PT = []
    for i in range(nb_blocks-6):
        state = F(state, 4)
        sr = get_sr(state)
        for j in range (64):
            recovered_PT.append(int(sr[j]) ^ int(CT_blocks[6+i][j]))
            sr[j] = CT_blocks[6+i][j]
        state = set_sr(state, sr)

    return recovered_PT


# **************************************
# **************************************
# Main function
# **************************************
# **************************************
def main():
    key = hexa_to_binary("2FAE237AEACBCD71AD5CC45D9F3649E4")
    IV  = hexa_to_binary("4CE5271616A805D6085E5F9E7BF2ECD3")
    # No associated data needed for the attacks
    AD  = ""
    PT  = hexa_to_binary("D201A35EE76A36A1CE847F5A9854E3ABD32B8F5546B4B116DEBB38B38A0977C86AE76C8B27BE05E2FADEA3F85DCBF23366CF75F29CDFAF6850AAFC2E87D214EB")
    hashlen = 256
    CT = hexa_to_binary("26abf3b5fcb254b4a7c0cae37f956130cf6e429a2bfedf6e6e3e9567854054d5f2bc1e2a5a7186758cccba24aea209138246df7c6fb9bc736a03a2bc894068fa")
    tag = hexa_to_binary("bdbebeac8992ef2d55866fc463b6fbe0233cf93decb72e7ecaa9d52320594de7")
    print(binary_to_hexa(Key_recovery_attack(PT, CT)))

    print("*****************************************")
    print("**** KEY ****")
    print('  '.join(elt for elt in binary_to_hexa(key).upper()))
    print("*****************************************")

    print("*****************************************")
    print("**** IV ****")
    print('  '.join(elt for elt in binary_to_hexa(IV).upper()))
    print("*****************************************")

    print("*****************************************")
    print("**** PLAINTEXT ****")
    print('  '.join(elt for elt in binary_to_hexa(PT).upper()))
    print("*****************************************")

    CT, tag = Encryption(key, IV, AD, PT, hashlen)
    print("*****************************************")
    print("ENCRYPTION DONE! WE GET THE FOLLOWING:")
    print("**** CIPHERTEXT ****")
    print('  '.join(elt for elt in binary_to_hexa(CT).upper()))
    print("*****************************************")

    print("*****************************************")
    print("**** TAG ****")
    print('  '.join(elt for elt in binary_to_hexa(tag).upper()))
    print("*****************************************")

    CT = ''.join(str(elt) for elt in CT)
    tag = ''.join(str(elt) for elt in tag)

    # Attacks
    recovered_key = Key_recovery_attack(PT, CT)
    recovered_pt = Plaintext_recovery_attack(CT)
    print("*****************************************")
    print("**** KEY RECOVERY ****")
    print('  '.join(elt for elt in binary_to_hexa(recovered_key).upper()))
    print("*****************************************")

    print("*****************************************")
    print("**** PLAINTEXT RECOVERY ****")
    print('  '.join(elt for elt in binary_to_hexa(recovered_pt).upper()))
    print("*****************************************")


main()