# ************************
# AUTHORS: Serkan MUTLU, Banu GUNER.
# Department of Industrial Engineering, Eskisehir Technical University, Eskisehir, Turkey.
# ************************
# Mixed-Model Two-Sided Disassembly Line Balancing (MTDLB) problem optimization codes by GUROBI.
# ************************
# 28th CIRP Life Cycle Engineering Conference - Jaipur, India.
# ************************


# -------------------------
# Import Library and Python Functions
# -------------------------


import gurobipy as grb
import time
import inputs


# -------------------------
# Inputs
# -------------------------


#PRE(m, k, i):      immediate normal node i predecessors of artificial node k for each model m
#SUC(m, k, i):      immediate normal node i successors of artificial node k for each model m
#THETA(m, i, s):    station-sides s where normal node i can be made for each model m
#C(m):              cycle time for each model m
#t(m,i):            processing time of each normal task i for each model m
print("MODELS\nModel-01\tFlashlight\nModel-02\tRadio\nModel-03\tToy Car\nModel-04\tBall Point Pen\n")
nA, nN, MODELS, STATIONS, SIDES, MONO, MONONO, PRE, SUC, THETA, C, t = inputs.MTDLBInput(int(input("If model-01 is disassemble in line -- 1; otherwise -- 0 = ")), int(input("If model-02 is disassemble in line -- 1; otherwise -- 0 = ")), int(input("If model-03 is disassemble in line -- 1; otherwise -- 0 = ")), int(input("If model-04 is disassemble in line -- 1; otherwise -- 0 = ")))

#EPSILON:           a very small number
#BIG_M:             a very big number
EPSILON = 0.001
BIG_M = 99999


# -------------------------
# Start of Gurobi Model
# -------------------------


time1 = time.time()
opt_model = grb.Model(name="MILP Model")


# -------------------------
# Variables -- MTDLB Variables
# -------------------------


# F_j    -- if station j is opened from both sides, 1; otherwise, 0.
F = opt_model.addVars(STATIONS, vtype=grb.GRB.BINARY)

# G_j    -- if station j is opened from only one sides, 1; otherwise, 0.
G = opt_model.addVars(STATIONS, vtype=grb.GRB.BINARY)

# U_js   -- if station j is opened from side s, 1; otherwise, 0.
U = opt_model.addVars(STATIONS, SIDES, vtype=grb.GRB.BINARY)

# Z_mi   -- if normal node i of model m is selected, 1; otherwise, 0.
Z = opt_model.addVars(MONO, vtype=grb.GRB.BINARY)


GAMMA = opt_model.addVars(MODELS, STATIONS, SIDES, vtype=grb.GRB.BINARY)

X = opt_model.addVars(MONO, STATIONS, SIDES, vtype=grb.GRB.BINARY)

DELTA = opt_model.addVars(MONONO, vtype=grb.GRB.BINARY)

TF = opt_model.addVars(MONO, lb=0)

# -------------------------
# Objective Function
# -------------------------

objective = (grb.quicksum(F[j] + G[j] for j in STATIONS) 
             + EPSILON*grb.quicksum(j*U[j,s] for j in STATIONS for s in SIDES)
             + EPSILON*EPSILON*grb.quicksum(t[m,i]*Z[m,i] for (m,i) in MONO))

# -------------------------
# Constraints
# -------------------------

# Equation (2)
opt_model.addConstrs(grb.quicksum(Z[m,i] for i in SUC[m,0]) == 1 for m in MODELS)

# Equation (3)
opt_model.addConstrs(grb.quicksum(Z[m,i] for i in SUC[m,k]) == grb.quicksum(Z[m,i] for i in PRE[m,k]) for m in MODELS for k in range(1, nA[m]) if ((m,k) in SUC) and ((m,k) in PRE))

# Equation (4)
opt_model.addConstrs(grb.quicksum(grb.quicksum(X[m,i,j,s] for s in THETA[m,i]) for j in STATIONS) == Z[m,i] for (m,i) in MONO)

# Equation (5)
opt_model.addConstrs(grb.quicksum(X[m,i,v,s] for i in PRE[m,k] for s in THETA[m,i] for v in range(1, j+1)) >= grb.quicksum(X[m,i,j,s] for i in SUC[m,k] for s in THETA[m,i]) for m in MODELS for k in range(1, nA[m]) for j in STATIONS if (((m,k) in PRE) and ((m,k) in SUC)))

# Equation (6)
opt_model.addConstrs(TF[m,i] <= C*Z[m,i] for (m,i) in MONO)

# Equation (7)
opt_model.addConstrs(TF[m,i] >= t[m,i]*Z[m,i] for (m,i) in MONO)

# Equation (8)
opt_model.addConstrs(TF[m,i] - TF[m,h] + BIG_M*(2 - grb.quicksum(X[m,i,j,s] for s in THETA[m,i]) - grb.quicksum(X[m,h,j,s] for s in THETA[m,h])) >= t[m,i] for m in MODELS for j in STATIONS for k in range(1, nA[m]) if (((m,k) in SUC) and ((m,k) in PRE)) for i in SUC[m,k] for h in PRE[m,k])

# Equation (9)
opt_model.addConstrs(TF[m,i] - TF[m,h] + BIG_M*(3 - X[m,i,j,s] - X[m,h,j,s] - DELTA[m,i,h]) >= t[m,i] for (m,i,h) in MONONO for j in STATIONS for s in SIDES)

# Equation (10)
opt_model.addConstrs(TF[m,h] - TF[m,i] + BIG_M*(2 - X[m,i,j,s] - X[m,h,j,s] + DELTA[m,i,h]) >= t[m,h] for (m,i,h) in MONONO for j in STATIONS for s in SIDES)

# Equation (11)
opt_model.addConstrs(grb.quicksum(X[m,i,j,s] for i in range(1, nN[m]+1) if s in THETA[m,i]) - nN[m]*GAMMA[m,j,s] <= 0 for m in MODELS for j in STATIONS for s in SIDES)

# Equation (12)
opt_model.addConstrs(grb.quicksum(GAMMA[m,j,s] for m in MODELS) - len(MODELS)*U[j,s] <= 0 for j in STATIONS for s in SIDES)

# Equation (13)
opt_model.addConstrs(grb.quicksum(U[j,s] for s in SIDES) - 2*F[j] - G[j] == 0 for j in STATIONS)

#-------------------------
# Model Setup
# -------------------------

opt_model.ModelSense = grb.GRB.MINIMIZE
opt_model.setObjective(objective)
opt_model.update()
opt_model.optimize()

print(opt_model)

print("\nSolution Results\n")
print("Time = ", time.time() - time1, "second")
print("Total number of stations opened from both sides\t\t:\t", sum([F[j].X for j in STATIONS]))
print("Total number of stations opened from only one side\t:\t", sum([G[j].X for j in STATIONS]))
print("Total number of stations opened\t\t\t\t:\t", sum([U[j,s].X for j in STATIONS for s in SIDES]))
for m in MODELS:
  print("#### MODEL-",m," ####")
  print("(m, i)\t\t (j,s)\t\t Processing Time\t Starting Time\t Ending Time")
  for i in range(1, nN[m]+1):
    if TF[m,i].X != 0.0:
      if i < 10:
        print((m,i), " :\t", [(j,s) for j in STATIONS for s in SIDES if X[m,i,j,s].X == 1.00], "\t", t[m,i], "\t\t\t", TF[m,i].X - t[m,i], "\t\t", TF[m,i].X)
      else:
        print((m,i), ":\t", [(j,s) for j in STATIONS for s in SIDES if X[m,i,j,s].X == 1.00], "\t", t[m,i], "\t\t\t", TF[m,i].X - t[m,i], "\t\t", TF[m,i].X)