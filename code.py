import math
import numpy as np
import matplotlib.pyplot as plt

def initialize():
    global L
    global t
    global n_of_webs
    global n_plates_top
    global n_plates_bottom
    global b_top
    global b_bot
    global first_car_weight
    global second_car_weight
    global locomotive_weight
    global dist_wheels_car
    global dist_cars_train
    global depth
    global d_web_plates
    global support_reaction_A
    global support_reaction_B
    global I_manual
    global y_centroid_manual
    global b_support
    global E
    global mu
    global b_in
    global b_out
    global a
    global loads
    global positions
    global tensile_strength
    global compressive_strength
    global shear_strength
    global shear_strength_cement
    global I_current
    global y_current
    global x_start
    global glue_tabs

    L = 1250 #mm
    t = 1.27 #mm
    n_of_webs = 2
    n_plates_top = 2
    n_plates_bottom = 1
    b_top = 100 #mm
    b_bot = 80 #mm
    depth = 100 #mm
    first_car_weight = 135 #N
    second_car_weight = 135 #N
    locomotive_weight = 182 #N
    dist_wheels_car = 176 #mm
    dist_cars_train = 164 #mm
    d_web_plates = 25 #mm
    support_reaction_A = 0
    support_reaction_B = 0
    I_manual = 85e3 #mm^4
    y_centroid_manual = 42.51 #mm
    b_support = 50 #mm
    E = 4000 #MPa N/mm^2
    mu = 0.2
    b_in = 100 #mm
    b_out = 0 #mm
    a=100 #mm
    tensile_strength = 30 #MPa
    compressive_strength = 6 #MPa
    shear_strength = 4 #MPa
    shear_strength_cement = 2 #MPa
    I_current = 0
    y_current = 0
    x_start = 172 #mm Center at 172 mm, min at 0 mm, max at L-908 mm
    glue_tabs = 10 #mm

def shear_strength_cement_calc():
    global shear_strength_cement
    Q = n_plates_top * b_top * t * (depth - n_plates_top*t/2 - y_current)
    #print("Q in shear strength of cement calc =", Q)
    shear_at_glue = (max_sfd*Q)/(n_of_webs*glue_tabs*I_current)    
    return shear_at_glue
    #print("Shear strength of cement set to:", shear_strength_cement, "MPa")


def max_bending_all_x():
    max_bending = 0.0
    for x in range(0, L-908 +1, 1):
        x_start = x
        train_loads(x_start)
        support_reactions()
        if max_BMD(L, loads, positions, support_reaction_A, support_reaction_B) > max_bending:
            max_bending = max_BMD(L, loads, positions, support_reaction_A, support_reaction_B)
    return max_bending
def max_shear_all_x():
    max_shear = 0.0
    for x in range(0, L-908 +1, 1):
        x_start = x
        train_loads(x_start)
        support_reactions()
        SFD_diagram(L, loads, positions, support_reaction_A, support_reaction_B)
        if max_sfd > max_shear:
            max_shear = max_sfd
    return max_shear


def desired_weight(N):
    global first_car_weight
    global second_car_weight
    global locomotive_weight
    second_car_weight = N/(3.618)
    first_car_weight = second_car_weight * 1.1
    locomotive_weight = first_car_weight * 1.38
    #print("First car weight set to:", first_car_weight, "N")
    #print("Second car weight set to:", second_car_weight, "N")
    #print("Locomotive weight set to:", locomotive_weight, "N")


def type_1():
    global E
    global mu
    global t
    global b_in
    stress = 0
    if b_in != 0:
        stress = ((4*math.pi**2*E)/(12*(1-mu**2)))*(t*n_plates_top/b_in)**2
    return stress

def type_2():
    global E
    global mu
    global t
    global b_out
    stress = 0
    if b_out != 0:
        stress = ((0.425*math.pi**2*E)/(12*(1-mu**2)))*(t*n_plates_top/b_out)**2
    return stress

def type_3():
    global E
    global mu
    global t
    distance = depth - n_plates_top*t - y_current
    stress = ((6*math.pi**2*E)/(12*(1-mu**2)))*(t/distance)**2
    return stress

def type_4():
    global E
    global mu
    global t
    global b_in
    global a
    stress = ((5*math.pi**2*E)/(12*(1-mu**2)))*((t/(depth-n_plates_top*t-n_plates_bottom*t))**2+(t/a)**2)
    return stress

def flexural_stress(I):
    
    # ensure loads/positions/supp reactions are set
    global n_of_webs, n_plates_top, n_plates_bottom, t, b_top, b_bot, depth
    # get centroid and section I if caller didn't pass an I

    # get signed bending moment with largest magnitude from max_BMD
    max_abs_M = max_BMD(L, loads, positions, support_reaction_A, support_reaction_B)
    # keep sign of the largest magnitude moment
    M =  max_abs_M

    # distances to extreme fibers from centroid (mm)
    y_bot = float(y_current)
    y_top = float(depth) - float(y_current)

    # bending stress = M * y / I  (N/mm^2 = MPa)
    sigma_top = M * y_top / I
    sigma_bot = M * y_bot / I

    # return (sigma_top_MPa, sigma_bot_MPa)
    return sigma_top, sigma_bot

def Q_calculation():
    A_top = b_top * n_plates_top*t
    A_web =  n_of_webs*t * (depth - t*(n_plates_top) - y_current)
    d_web = (depth - t*n_plates_top - y_current)/2
    d_top = (depth - n_plates_top*t/2 - y_current) 
    print(depth, t, n_plates_top, y_current)
    print("d_web =", d_web, "d_top =", d_top, "A_web =", A_web, "A_top =", A_top)
    Q = A_top * d_top + A_web * d_web 
    print("Q =", Q)
    return Q

def x_for_mid_location(L):
    train_len = 3*dist_wheels_car + 2*dist_cars_train
    return (L-train_len)/2 - b_support/2

def shear_stress(I):
    
    # ensure loads/positions/supp reactions are set
    global n_of_webs, n_plates_top, n_plates_bottom, t, b_top, b_bot, depth
    # get centroid and section I if caller didn't pass an I

    # get shear force with largest magnitude from SFD
    max_shear = max_sfd

    # Q calculation
    Q = Q_calculation()

    # shear stress = V * Q / (I * t)  (N/mm^2 = MPa)
    tau_max = max_shear * Q / (I * t*n_of_webs)

    return tau_max

def train_loads(x_start):
    global loads, positions, first_car_weight, second_car_weight, locomotive_weight
    #print("Setting train loads with first car weight:", first_car_weight, "N, second car weight:", second_car_weight, "N, locomotive weight:", locomotive_weight, "N")
    loads = [first_car_weight/2, first_car_weight/2, second_car_weight/2, second_car_weight/2, locomotive_weight/2, locomotive_weight/2]
    positions = [x_start, x_start + dist_wheels_car, x_start + dist_wheels_car + dist_cars_train, x_start + 2*dist_wheels_car + dist_cars_train, x_start + 2*dist_wheels_car + 2*dist_cars_train, x_start + 3*dist_wheels_car + 2*dist_cars_train]
    # Locomotive
    #print("Location of loads (mm):", positions)


def support_reactions():
    global support_reaction_A
    global support_reaction_B
    total_moment = 0.0
    for i in range(len(loads)):
        total_moment += loads[i] * (positions[i])
    support_reaction_B = total_moment / 1200
    #print("Total moment about A:", total_moment, "N*mm", "Support reaction B calculated as:", support_reaction_B, "N")
    #print("Total loads:", sum(loads), "N")
    support_reaction_A = sum(loads) - support_reaction_B


def I_calculation(n_of_webs, n_plates_top, n_plates_bottom, t, b_top, b_bot, depth):
    
    # web free height (remove both top and bottom plate thicknesses)
    h_web = float(depth) - (n_plates_top + n_plates_bottom) * float(t)
    #print("h_web =", h_web)
    if h_web <= 0:
        raise ValueError("Invalid geometry: web height <= 0")

    # areas (mm^2)
    A_web = float(n_of_webs) * float(t) * h_web
    A_top = float(n_plates_top) * float(b_top) * float(t)
    A_bot = float(n_plates_bottom) * float(b_bot) * float(t)
    #print("A_web =", A_web, "A_top =", A_top, "A_bot =", A_bot)
    A_total = A_web + A_top + A_bot
    if A_total == 0.0:
        raise ValueError("Total area is zero")

    # y positions from bottom datum (mm)
    y_bot = n_plates_bottom*float(t) / 2.0
    # web centroid sits above the bottom plates
    y_web = (float(n_plates_bottom) * float(t) + h_web/2 )
    y_top = depth - (n_plates_top*float(t) / 2.0)
    #print("y_bot =", y_bot, "y_web =", y_web, "y_top =", y_top)
    # centroid (mm)
    y_centroid = (A_bot * y_bot + A_web * y_web + A_top * y_top) / A_total
    #print("y_centroid =", y_centroid)
    # individual second moments about their own centroids (mm^4)
    I_web_ind = (n_of_webs*float(t) * h_web**3) / 12.0
    I_top_ind = (float(b_top) * (n_plates_top*float(t))**3) / 12.0
    I_bot_ind = (float(b_bot) * (n_plates_bottom*float(t))**3) / 12.0
    #print("I_web_ind =", I_web_ind, "I_top_ind =", I_top_ind, "I_bot_ind =", I_bot_ind)
    # total I using parallel axis theorem
    #print(y_top - y_centroid)
    I_web = (I_web_ind + (A_web* (y_centroid - y_web)**2))
    I_top =  (I_top_ind + ((A_top) * (y_top - y_centroid)**2))
    I_bot =  (I_bot_ind + (A_bot) * (y_centroid - y_bot)**2)
    #print("I_web =", I_web, "I_top =", I_top, "I_bot =", I_bot)



    I_total = I_web + I_top + I_bot

    return y_centroid, I_total



def SFD_diagram(L, loads, positions, support_reaction_A, support_reaction_B):
    global max_sfd

    shear = support_reaction_A
    max_sfd = abs(shear)

    for load in loads:
        shear = shear - load
        if abs(shear) > max_sfd:
            max_sfd = abs(shear)
    
    shear = shear + support_reaction_B
    if abs(shear) > max_sfd:
        max_sfd = abs(shear)
    return max_sfd


def max_BMD(L, loads, positions, support_reaction_A, support_reaction_B):
    global max_bending_moment

    shear = support_reaction_A
    moment = 0
    max_bending_moment = 0

    for i in range(len(loads)):
        load = loads[i]
        pos = positions[i]
        dist = pos - (0 if i == 0 else positions[i-1])
        moment += shear * dist
        if abs(moment) > abs(max_bending_moment):
            max_bending_moment = moment
        shear -= load


    last_dist = L - positions[-1]
    moment += shear * last_dist
    if abs(moment) > abs(max_bending_moment):
        max_bending_moment = moment
    
    return abs(max_bending_moment)




    # L = float(L)
    # # event points to evaluate exactly (supports and load positions)
    # events = [0.0, L]
    # events += [float(p) for p in positions]
    # events = sorted(set(min(max(e, 0.0), L) for e in events))

    # def M_at(x):
    #     # bending moment about left support (x=0)
    #     M = float(support_reaction_A) * x
    #     for p, w in zip(positions, loads):
    #         if x > p:
    #             M -= w * (x - p)
    #     return M

    # # sample uniform grid plus event points
    # x_samples = [L * i / samples for i in range(samples + 1)]
    # for e in events:
    #     if e not in x_samples:
    #         x_samples.append(e)
    # x_samples = sorted(x_samples)

    # M_samples = [M_at(x) for x in x_samples]
    # max_pos = max(M_samples)
    # max_neg = min(M_samples)
    # abs_vals = [abs(m) for m in M_samples]
    # max_abs = max(abs_vals)
    # idx = abs_vals.index(max_abs)
    # x_at_max = x_samples[idx]
    # #print("Max Bending Moment:", max_abs, "at x =", x_at_max)
    # return max_abs, x_at_max, max_pos, max_neg
# Function that will plot all the required figures
def sfe():
    global sfe_shears
    global locations_sfe
    sfe_shears = []
    locations_sfe = []
    bme_moments = []
    factor_safety = []
    for x in range(0, L-908 +1, 1):
        x_start = x
        train_loads(x_start)
        support_reactions()
        sfe_shears.append(SFD_diagram(L, loads, positions, support_reaction_A, support_reaction_B))
        bme_moments.append(max_BMD(L, loads, positions, support_reaction_A, support_reaction_B))
        locations_sfe.append(x)
        factor_safety.append(check_safety())
    plt.figure(1)
    plt.subplot(211)
    plt.plot(locations_sfe, sfe_shears)
    plt.ylabel("Shear Force (N)")
    plt.title("Shear Force Envelope along Bridge Length")
    plt.xlabel("Location along bridge (mm)")
    plt.subplot(212)
    plt.plot(locations_sfe, bme_moments)
    plt.xlabel("Location along bridge (mm)")
    plt.ylabel("Bending Moment (N*mm)")
    plt.title("Bending Moment Envelope along Bridge Length")
    plt.grid()
    plt.show()

    plt.figure(2)
    plt.plot(locations_sfe, factor_safety)
    plt.xlabel("Location along bridge (mm)")
    plt.ylabel("Load (N)")
    plt.title("Max load before failure along Bridge Length")
    plt.grid()
    plt.show()

        


def print_geometry():
    print("Geometry:")
    print(f"  n_of_webs = {n_of_webs}")
    print(f"  n_plates_top = {n_plates_top}")
    print(f"  n_plates_bottom = {n_plates_bottom}")
    print(f"  t = {t} mm")
    print(f"  b_top = {b_top} mm")
    print(f"  b_bot = {b_bot} mm")
    print(f"  depth = {depth} mm")

def results():
    train_loads(x_start)
    support_reactions()
    print("Support A:", support_reaction_A, "Support B:", support_reaction_B)
    SFD_diagram(L, loads, positions, support_reaction_A, support_reaction_B)
    print("Max Shear Force (N):", max_sfd)
    max_BMD(L, loads, positions, support_reaction_A, support_reaction_B)
    print("Max Bending Moment (N*mm):", max_BMD(L, loads, positions, support_reaction_A, support_reaction_B))
    print("Flexural stresses (MPa) at top and bottom:", flexural_stress(I_current))
    #print("Q =", Q_calculation())
    print("Shear stress:", shear_stress(I_current)) 
    print("Shear strength of cement (MPa):", shear_strength_cement_calc())
    print("Type 1 buckling stress (MPa):", type_1())
    print("Type 2 buckling stress (MPa):", type_2())
    print("Type 3 buckling stress (MPa):", type_3())
    print("Type 4 buckling stress (MPa):", type_4()) 

def check_safety():
    sigma_top, sigma_bot = flexural_stress(I_current)
    tau_max = shear_stress(I_current)
    fail_factor_of_safety = 100
    safe = True
    if abs(sigma_top) > compressive_strength:
        print(f"\033[91mFail: Top flexural stress {sigma_top:.2f} MPa exceeds compressive strength {compressive_strength} MPa\033[0m")
        safe = False
        if compressive_strength/abs(sigma_top) < fail_factor_of_safety:
            fail_factor_of_safety = compressive_strength/abs(sigma_top)
    else:
        print("Safety facrtor for compression:", compressive_strength/abs(sigma_top), "\033[92mSafe \033[91m\033[0m")
    if abs(sigma_bot) > tensile_strength:
        print(f"\033[91mFail: Bottom flexural stress {sigma_bot:.2f} MPa exceeds tensile strength {tensile_strength} MPa\033[0m")
        safe = False
        if tensile_strength/abs(sigma_bot) < fail_factor_of_safety:
            fail_factor_of_safety = tensile_strength/abs(sigma_bot)
    else:
        print("Safety facrtor for tension:", tensile_strength/abs(sigma_bot), "\033[92mSafe \033[91m\033[0m")
    if abs(tau_max) > shear_strength:
        print(f"\033[91mFail: Shear stress {tau_max:.2f} MPa exceeds shear strength {shear_strength} MPa\033[0m")
        safe = False
        if shear_strength/abs(tau_max) < fail_factor_of_safety:
            fail_factor_of_safety = shear_strength/abs(tau_max)
    else:
        print("Safety facrtor for shear:", shear_strength/abs(tau_max), "\033[92mSafe \033[91m\033[0m")
    if shear_strength_cement_calc() > shear_strength_cement:
        print(f"\033[91mFail: Shear at glue {shear_strength_cement_calc():.2f} MPa exceeds shear strength of cement {shear_strength_cement} MPa\033[0m")
        safe = False
        if shear_strength_cement_calc()/shear_strength_cement < fail_factor_of_safety:
            fail_factor_of_safety = shear_strength_cement_calc()/shear_strength_cement
    else:
        print("Safety facrtor for shear at glue:", shear_strength_cement/shear_strength_cement_calc(), "\033[92mSafe \033[91m\033[0m")
    if type_1() < sigma_top and type_1() != 0:
        print(f"\033[91mFail: Type 1 buckling stress {type_1():.2f} MPa less than compressive strength {sigma_top} MPa\033[0m")
        safe = False
        if type_1()/sigma_top < fail_factor_of_safety:
            fail_factor_of_safety = type_1()/sigma_top
    else:
        print("Safety facrtor for Type 1 buckling:", type_1()/sigma_top, "\033[92mSafe \033[91m\033[0m")
    if type_2() < sigma_top and type_2() != 0:
        print(f"\033[91mFail: Type 2 buckling stress {type_2():.2f} MPa less than compressive strength {sigma_top} MPa\033[0m")
        safe = False
        if type_2()/sigma_top < fail_factor_of_safety:
            fail_factor_of_safety = type_2()/sigma_top
    else:
        print("Safety facrtor for Type 2 buckling:", type_2()/sigma_top, "\033[92mSafe \033[91m\033[0m")
    if type_3() < sigma_top:
        print(f"\033[91mFail: Type 3 buckling stress {type_3():.2f} MPa less than compressive strength {sigma_top} MPa\033[0m")
        safe = False
        if type_3()/sigma_top < fail_factor_of_safety:
            fail_factor_of_safety = type_3()/sigma_top
    else:
        print("Safety facrtor for Type 3 buckling:", type_3()/sigma_top, "\033[92mSafe \033[91m\033[0m")
    if type_4() < tau_max:
        print(f"\033[91mFail: Type 4 buckling stress {type_4():.2f} MPa less than compressive strength {tau_max} MPa\033[0m")
        safe = False
        if type_4()/tau_max < fail_factor_of_safety:
            fail_factor_of_safety = type_4()/tau_max
    else:
        print("Safety facrtor for Type 4 buckling:", type_4()/tau_max, "\033[92mSafe \033[91m\033[0m")
    if safe:
        print("\033[92mAll checks passed: Structure is safe.\033[0m")
    if not safe:
        print(f"\033[91mSome checks failed: Structure is NOT safe. Lowest factor of safety: {fail_factor_of_safety:.2f}, it will fail if loads are increased by a weight of {fail_factor_of_safety * 2000:.2f} N or more.\033[0m")
    return fail_factor_of_safety*2000
print("\033[92mSafe \033[91mFail\033[0m")
#print("First car weight:", first_car_weight, "N", "Second car weight:", second_car_weight, "N", "Locomotive weight:", locomotive_weight, "N")
print("\033[93m ------------ TESTING FOR BASE LOAD CASE 2 ------------\033[00m")
initialize()
print("Testing at location x_start =", x_start, "mm")
desired_weight(2000)
I_current = I_calculation(n_of_webs, n_plates_top, n_plates_bottom, t, b_top, b_bot, depth)[1]
y_current = I_calculation(n_of_webs, n_plates_top, n_plates_bottom, t, b_top, b_bot, depth)[0]
print("I_current =", I_current, "mm^4", "y_current =", y_current, "mm")
results()
print("\033[93m ------------ SAFETY CHECKS -------------\033[00m")
check_safety()
print("Max bending moment at all x (N*mm):", max_bending_all_x())
print("Max shear force at all x (N):", max_shear_all_x())
sfe()
#print(I_calculation(n_of_webs, n_plates_top, n_plates_bottom, t, b_top, b_bot, depth))
#print(I_manual, y_centroid_manual)
# y, I = I_calculation(n_of_webs, n_plates_top, n_plates_bottom, t, b_top, b_bot, depth)
# print(f"y_centroid = {y:.6f} mm")
# print(f"I_total  = {I:.6f} mm^4  ({I/1e4:.6f} x10^4 mm^4)")
# train_loads(172)
# support_reactions()
# print(x_for_mid_location(L))
# print("Support A:", support_reaction_A, "Support B:", support_reaction_B)
# SFD_diagram(L, loads, positions, support_reaction_A, support_reaction_B)
# max_BMD(L, loads, positions, support_reaction_A, support_reaction_B)
# print("Flexural stresses (MPa) at top and bottom:", flexural_stress(I_manual))
# print("Q =", Q_calculation())
# print("Shear stress:", shear_stress(I_manual))