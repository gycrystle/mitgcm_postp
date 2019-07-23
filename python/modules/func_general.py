def profile_gen(x, a, b, c, depth, gamma):
  return a * np.exp(-b * (x/depth)**gamma) + c
