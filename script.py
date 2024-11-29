import os
st = [2*i for i in range(3, 26)]
# sparsity_approx = []
# time_approx = []
# sparsity_exact = []
# time_exact = []
# for n in st:
#     os.system(f"python3 graphgen.py {n}")
#     s = os.popen(f"python3 frechet_code.py").read()
#     sparsity_approx.append(round(float(s.split()[0]), 2))
#     time_approx.append(round(float(s.split()[1]), 2))
#     t = os.popen(f"python3 scip_demands.py").read()
#     sparsity_exact.append(round(float(t.split()[0]), 2))
#     time_exact.append(round(float(t.split()[1]), 2))
#     print(sparsity_approx)
#     print(time_approx)
#     print(sparsity_exact)
#     print(time_exact)
#     print()

# print("_________________________________________")
# print(sparsity_approx)
# print(time_approx)
# print(sparsity_exact)
# print(time_exact)



# [0.7, 0.48, 0.49, 0.22, 0.27, 0.36, 0.33, 0.39, 0.38, 0.37]
# [0.02, 0.08, 0.29, 0.86, 1.99, 3.5, 6.49, 10.04, 14.78, 22.05]
# [0.7, 0.48, 0.48, 0.22, 0.27, 0.35, 0.35, 0.33, 0.35, 0.37]
# [0.26, 0.52, 0.81, 0.53, 0.54, 7.12, 12.09, 36.8, 158.99, 271.38]


sparsity_approx =  [0.36, 0.42, 0.48, 0.48, 0.42, 0.56, 0.37, 0.22, 0.28, 0.28, 0.25, 0.4, 0.36, 0.35, 0.38, 0.38, 0.32, 0.37, 0.25, 0.3, 0.33, 0.31, 0.37]
time_approx = [0.03, 0.05, 0.13, 0.21, 0.39, 0.79, 1.03, 1.37, 1.82, 2.34, 3.57, 4.58, 5.46, 7.19, 8.22, 11.65, 13.2, 15.89, 18.37, 25.43, 24.59, 30.15, 23.27]
sparsity_exact = [0.36, 0.42, 0.48, 0.48, 0.42, 0.48, 0.37, 0.22, 0.28, 0.27, 0.25, 0.4, 0.35, 0.35, 0.38, 0.37, 0.32, 0.33, 0.25, 0.3, 0.33, 0.31, 0.37]
time_exact = [0.35, 0.24, 0.74, 0.59, 1.01, 0.87, 1.32, 1.05, 1.66, 1.48, 3.66, 6.46, 11.41, 9.63, 17.79, 81.98, 26.75, 59.16, 46.63, 70.96, 179.25, 91.96, 560.02]
sparsity_ratio = [sparsity_approx[i]/sparsity_exact[i] for i in range(len(sparsity_approx))]
import matplotlib.pyplot as plt
import math

# # Plot of fraction of sparsity vs number of nodes, along with a plot of log(n)
# plt.figure(figsize=(10, 6))
# plt.plot(st, sparsity_ratio, label="Sparsity ratio")
# # log n
# # plt.plot(st, [1+0.1*(math.log(n)-1) for n in st], label="log(n)")
# plt.xlabel("Number of nodes")
# plt.ylabel("Sparsity ratio")
# plt.title("Sparsity ratio vs Number of nodes")
# plt.legend()
# plt.grid()
# plt.show()

# Plot the time taken for the approximate and exact methods
plt.figure(figsize=(10, 6))
plt.plot(st, time_approx, label="Approximate method")
plt.plot(st, time_exact, label="Exact method")
plt.xlabel("Number of nodes")
plt.ylabel("Time taken (s)")
plt.title("Time taken for approximate and exact methods")
plt.legend()
plt.grid()
plt.show()


