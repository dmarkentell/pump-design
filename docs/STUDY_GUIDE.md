# Study guide: the reasoning behind every decision

This is the long version of the README. For each step it explains what was done, **why**, what
the equation means physically, and what would change if the choice were different. The last
section is a self-quiz with the questions an interviewer is most likely to ask.

Numbers are for the example duty (660 gpm, 295 ft, 3560 rpm) and the optimized impeller unless
stated. File references point to where each idea lives in the code.

---

## 0. The big picture: what depends on what

```
duty point (Q, H, n, NPSHa)
   │
   ├─► motor size ─► shaft torque ─► shaft Ø at impeller ─► hub Ø        (shaftrule.py)
   │                                                          │
   └─────────────────────────────────────────────────────────► impeller + volute sizing  (sizing.py)
                                                                   │
                         1D performance model (meanline.py) ◄──────┘
                                   │
                    calibrate one loss factor (calibrate.py)
                                   │
                  sweeps + optimizer under screening checks (run_study.py)
                                   │
                         optimized impeller geometry
                         ├─► radial + axial loads ─► shaft, bearings, key, critical speed (drive_end.py)
                         └─► 3D solid ─► stress under spin + pressure (impeller_fea.py)
```

The order matters. The hub has to be known before the eye can be sized, because the hub blocks
part of the eye. The impeller has to be known before the drive end, because the impeller makes the
loads. **This ordering was a mistake in the first version**: the hub was a placeholder, and fixing
it changed the hydraulic results (section 8, problem 6).

---

## 1. Pump vocabulary you must know cold

**Head (H).** Energy added per unit weight of fluid, in metres or feet. It is used instead of
pressure because a pump makes the same head whatever the fluid's density; pressure rise is
ρgH. 90 m of head is about 880 kPa for water.

**Flow (Q).** Volume per time. 150 m³/h = 0.0417 m³/s = 660 gpm.

**Best efficiency point (BEP).** The flow where efficiency peaks. Away from it, flow meets the
blades and volute at the wrong angles and the radial load grows. API 610 wants the rated flow
between **80% and 110% of BEP**, and calls **70–120%** the preferred operating region.

**Specific speed.** One number that describes the *shape* of the pump the duty calls for:
nq = n·√Q / H^0.75 (rpm, m³/s, m), which is 24.9 here. US units: Ns = n·√gpm / ft^0.75, which is 1285 here.
Low specific speed means high head and low flow: a large diameter and narrow, radial passages. That is
exactly this pump. Every statistical correlation in the model (head coefficient, outlet width,
efficiency) is a function of nq, because pumps with the same nq look alike.

**NPSH.** Net positive suction head: how far the suction pressure sits above the vapor pressure,
expressed as head. *Available* (NPSHa) is set by the system (10 m here). *Required* (NPSHr) is
set by the pump: the pressure drops as fluid accelerates into the eye, and if it drops to vapor
pressure, the fluid boils (cavitation). Rule: NPSHa must exceed NPSHr with margin (1 m here).

**Affinity laws.** For the same pump at a different speed: Q ∝ n, H ∝ n², P ∝ n³. For a trimmed
impeller at the same speed, roughly Q ∝ D, H ∝ D². The tests check the first set
(`test_affinity_laws`).

**OH2.** API 610's code for an *overhung* pump (impeller on the end of the shaft, beyond the
bearings), *centerline-mounted* (the casing feet are at the shaft centerline, so thermal growth
does not misalign it with hot fluids), with its own bearing housing and a flexible coupling to
the motor. It is the workhorse refinery pump.

---

## 2. The 1D hydraulic model (`pump/meanline.py`)

"1D" or "meanline" means the flow is represented by one average streamline through the impeller,
with velocities at the inlet and outlet. It cannot see 3D flow features, but it captures the energy
balance, which is what sets head, power and efficiency.

### 2.1 Velocity triangles

At any point on a blade there are three velocities:
- **u**, blade speed = ω·r. At the outlet, u2 = 41.6 m/s.
- **w**, velocity *relative to the blade* (what you would see riding on the impeller). It follows
  the blade surface.
- **c**, *absolute* velocity (what a stationary observer sees), c = u + w as vectors.

c splits into a **meridional** part cm (through-flow, set by Q and the flow area) and a
**tangential** part cu (swirl). At the outlet: cm2 = Q/(π·D2·b2) = 4.5 m/s.

### 2.2 Euler's equation: where the head comes from

H_th = (u2·cu2 − u1·cu1) / g. The impeller adds energy by changing the fluid's angular momentum.
With no swirl entering (cu1 = 0, an assumption for a plain suction with no inlet guide vanes),
H_th = u2·cu2 / g = 99.4 m at rated flow for the optimized impeller. **Everything else in the model
is about how much of that ideal head is lost.**

Why head falls as flow rises: more flow raises cm2, and with a backward-swept blade that lowers
cu2 = u2 − cm2/tan(β2). So even the ideal head line slopes down with flow. That slope is what makes
a stable curve possible. A blade angle of β2 = 90° (radial) would give a flat ideal line.

### 2.3 Slip (Wiesner correlation)

The fluid does not leave exactly along the blade. A rotating channel sets up a relative eddy, so
the flow lags the blade and cu2 is smaller than the blade angle implies. The slip factor
σ = 1 − √(sin β2) / z^0.7 (Wiesner), corrected when the blades are short relative to the impeller.
Here σ = 0.78: **slip costs about 22% of the ideal swirl**. More blades means less slip, which is
the one real advantage of more blades. Code: `_slip_wiesner`.

### 2.4 Blockage

Blades have thickness, so they take up flow area. The factor τ = 1/(1 − z·t/(π·D·sin β)) raises
the real through-flow velocity. It matters most at the eye, where it raises NPSHr. That is why
7+ blades fail the NPSH check.

### 2.5 The losses, and what each does to the curve

| Loss | Formula idea | Behavior | Why |
|---|---|---|---|
| Incidence | k·(Δwu)²/2g, where Δwu is the mismatch between flow angle and blade angle at the inlet | zero near design flow, grows both ways | flow hitting the blade at an angle separates. A softer coefficient below design (0.3 vs 0.7), because recirculation cushions it |
| Impeller friction | 4·cf·(L/Dh)·w²/2g, the pipe-friction form on the blade channel | grows as Q² | wall shear in the channel |
| Volute mismatch | k·(cu3 − c_throat)²/2g | zero where the volute matches the impeller | at off-design, the fluid leaving the impeller is too fast or too slow for the volute area and has to decelerate or accelerate |
| Volute friction | k·c_throat²/2g | grows as Q² | wall shear in the volute |
| Discharge cone | k·(c_throat − c_discharge)²/2g | grows as Q² | diffuser loss slowing the flow to the pipe |

At rated flow on the optimized impeller: friction 4.8 m, volute friction 2.4 m, cone 1.4 m,
mismatch 0.6 m, incidence 0.15 m. Hydraulic efficiency H/H_th = 90.5%.

**Why the BEP exists:** friction losses grow with Q², the shock-type losses (incidence, mismatch)
grow away from the design flow in both directions, and disk friction is a fixed tax that hurts
most at low flow. Efficiency peaks where these balance.

### 2.6 Leakage, and why it is iterated

Some flow goes round again: from the high-pressure side gap, through the wear-ring clearance,
back to the eye. That flow passes through the impeller twice, so it costs power but delivers
nothing. The driving pressure is the static pressure the impeller builds, minus the drop from the
fluid spinning in the side gap at about half impeller speed. **The circular dependency**: leakage
depends on impeller pressure, impeller pressure depends on the flow through the impeller, and that
flow includes the leakage. So the model iterates six times (it converges in 2–3). Volumetric
efficiency comes out at 94%, typical for API clearances (0.40 mm diametral), which are
deliberately generous so a rub on a shaft bow does not seize the pump.

### 2.7 Disk friction (Daily–Nece)

The outside of both shrouds spins in fluid, which takes power and adds no head. The moment
coefficient is Cm = 0.0255·(s/r)^0.1 / Re^0.2, power = Cm·ρ·ω³·r⁵. **It scales as D⁵**, which is
why low-specific-speed pumps (big diameter for their flow) lose efficiency to it. It is about 1.5%
of power here.

### 2.8 NPSH required, and why the eye has an optimum size

NPSHr = λc·cm1²/2g + λw·w1²/2g. The first term is the acceleration into the eye. The second is the
low-pressure peak on the blade's leading edge. A **bigger eye lowers cm1** (more area) but **raises
w1** (higher blade speed at a bigger radius). Minimizing the sum gives
D1 = 4.5·(Q / (k·n))^(1/3), where k = 1 − (dh/D1)² is the share of the eye the hub leaves free.
The sweep confirms the minimum near 1.0–1.05× this size. The optimizer goes *smaller* (0.93×)
because a smaller eye means a smaller wear ring, which means less leakage and better efficiency, until
NPSHr hits the margin.

### 2.9 How a design is completed (`pump/sizing.py`)

For any set of design choices, the code builds a *complete* pump so designs can be compared fairly:
1. **Eye and hub** from the rules above.
2. **Outlet diameter D2 is solved** (Brent's root finder) so the head at rated flow is exactly 90 m.
   This is the computer version of trimming an impeller to hit the rated point.
3. **Inlet blade angle is iterated** so the blade sits 3° steeper than the real flow angle at rated
   flow (+3° incidence, a standard choice: slight positive incidence protects against cavitation at
   higher flows). The loop is needed because blockage and leakage change the flow angle.
4. **Volute throat is sized** so the BEP lands at the rated flow. In real pumps the volute throat
   is the main control on where the BEP falls.

---

## 3. Calibration (`pump/calibrate.py`): why one factor, and why only friction

Textbook loss coefficients are generic. Rather than trust them, the model is anchored to real
data: Gülich's statistical correlations, fitted to hundreds of measured pumps, give the typical
hydraulic efficiency for a pump of this size and nq (89.7%). One multiplier on the friction losses
is solved so the baseline hits that number. It came out at **2.50**: the textbook friction values
alone were too optimistic.

**Why only the friction losses?** The calibration point is the BEP. At BEP, the shock-type losses
(incidence, mismatch) are near zero, so the data says nothing about them. Scaling them would just
distort the off-design curve shape with no evidence.

**Why the factor is frozen afterwards.** If it were re-tuned for every design, every design would
match the statistics by construction and the comparison would mean nothing. Frozen, the
differences between designs come from their geometry.

**The independent check.** Gülich also gives an *overall* efficiency correlation, which includes
leakage, disk friction and mechanical losses. That one was not used to calibrate. The model predicts
81.8%; the correlation says 82.0%. Agreement on a quantity you didn't tune is real evidence.

---

## 4. Optimization (`run_study.py`)

**Variables:** blade count z (5–9), outlet angle β2 (14–36°), outlet width (0.6–1.4× the
correlation), eye diameter (0.85–1.25× the NPSH-optimal size), inlet incidence (0–6°).

**Objective:** efficiency at rated flow.

**Constraints and why each exists:**
- NPSH margin ≥ 1 m: cavitation erodes impellers and kills performance.
- Rated flow 80–110% of BEP: API 610. Off-BEP operation means high radial load and recirculation.
- Head rise to shutoff ≥ 10%: with a flat curve, a small pressure change swings the flow a lot,
  so control is poor and parallel pumps fight each other.
- Rising curve (no hump): with a hump, two flows give the same head, and the pump can hunt between them.
- w2/w1 ≥ 0.7 (de Haller): too much deceleration of the relative flow and it separates off the blade.
- Outlet width ≥ 8 mm: a narrower passage is hard to cast cleanly.
- Blade count 5–9: see the Coppage story below.

**Method:** differential evolution (a population-based global search) for each blade count,
with constraint violations added as a penalty. A global search was chosen because the design space
has several local optima (the "ridge" below) and the constraint boundaries make the objective
non-smooth, so a gradient method would get stuck.

**Results and their meaning:**
- **+1.0 point** (81.8 → 82.8%), limited by NPSH margin. The optimizer always runs into the
  efficiency-vs-suction trade.
- **The optimum is a flat ridge.** A narrow outlet with a steep angle and a wider outlet with a
  shallower angle make the same head (head depends on cu2, which depends on both b2 and β2). The
  best 5- and 6-blade designs are within 0.2 points. That tells you the 1D model cannot pick
  between them; castability, pulsation and curve shape would.
- **Steep β2 makes the curve hump** from about 24° up. This is real pump behavior: a steep blade
  flattens the ideal head line (section 2.2), and the losses at low flow then bend it over.

**The Coppage experiment (a good interview story).** The optimizer first picked 4 blades. Fewer
blades means less friction and blockage, and the 1D model can't see why too few blades is bad
(poor flow guidance, high blade loading, strong pressure pulses). I added a blade-loading loss from
compressor practice (Coppage's diffusion factor). It made things worse: it put a 3% hump near
shutoff on a perfectly ordinary 6-blade design, because the correlation was built for the design
point. So I switched it off, kept the diffusion factor as a reported number, and bounded blade
count by practice instead. **The lesson: a correlation outside its range can do more harm than
leaving the effect out, and saying so is better than hiding it.**

---

## 5. Drive end (`pump/drive_end.py`)

### 5.1 Motor
API 610 requires the motor nameplate to be at least 125% / 115% / 110% of rated pump power, for
pumps up to 22 kW / 22–55 kW / above 55 kW. The margin covers wear, fluid changes and end-of-curve
operation. Rated power is 44 kW → 115% → 51 kW → the next NEMA size is **75 hp**. Check: power at
120% of BEP with the maximum impeller is 71.6 hp, still under 75.

### 5.2 Loads on the impeller

**Radial (Stepanoff).** A single volute is designed for one flow. Away from it, pressure around
the impeller is uneven and pushes it sideways: F = K·ρgH·D2·B2, with K = 0.36·|1 − (Q/Q_BEP)²|.
That is near zero at BEP (floored at 0.05, because nothing balances perfectly) and largest at
shutoff. **Why the max impeller and 40–120% BEP:** API 610 requires that the casing accept an impeller
giving 5% more head, so the shaft must survive that impeller, and the pump must survive its whole
allowable operating range, not just the rated point. Worst case is 1.81 kN; at rated flow it is only 0.24 kN.

**Axial (thrust).** The back shroud sees high pressure over more area than the front, because the eye
is at suction pressure. Balance holes plus a back wear ring at the same diameter as the front ring
equalize the two sides, leaving small terms: inlet momentum (pushes toward the motor), suction
pressure on the shaft's cross-section, and the small pressure drop across the balance holes. Real
balanced impellers still see a thrust nobody can predict well, because the fluid spins at
different speeds on the two sides. So I add ±10% of one shroud's force. Nominal thrust is 0.16 kN,
design thrust is 2.6 kN. The allowance dominates, and that is the honest answer: **you design a thrust
bearing for the uncertainty, not the nominal**.

### 5.3 Shaft sizing logic
1. **Impeller fit diameter** from torque: τ = 16T/(πd³) at 40 MPa design shear with the full motor
   torque (150 N·m) → 26.7 mm → **30 mm**. The low design stress leaves room for the keyway and fatigue.
2. **Seal section** starts at 40 mm and steps up 5 mm until **both** (a) deflection at the seal
   faces is ≤ 45 µm (50 µm limit with 10% margin) and (b) the inboard bearing shoulder passes
   fatigue. Result: **45 mm**.
3. Bearing journals are 5 mm larger (**50 mm**), and the span between bearings 10 mm more (**60 mm**),
   so each bearing can slide on over a smaller diameter and seat against a shoulder.

**Why 50 µm at the seal faces?** A mechanical seal's faces are lapped flat to within a few light bands.
If the shaft whips, the faces open, leak and wear. API 610 sets the limit at the seal, not at the
impeller, because the seal is what fails.

### 5.4 Beam finite-element model
The shaft is split into 4 mm Euler–Bernoulli beam elements, each with the local diameter. Each
node has a deflection and a slope. Bearings are springs. Solving K·u = F gives the deflected shape;
the same K with a mass matrix M gives the natural frequencies from K·φ = ω²·M·φ.
Checked against hand formulas: tip deflection of an overhung beam F·a²(a+L)/(3EI) to 0.1%, and the
natural frequency of a point mass on an overhang, √(k/m), to 0.2%.

**Bearing stiffness is the assumption that matters most.** Each bearing gives a little under load.
With soft bearings (≤1×10⁸ N/m), the bearings alone let the seal faces move more than 50 µm, and **no
shaft diameter fixes it**. At 3×10⁸ N/m (typical for preloaded ball bearings this size in a stiff
housing) fatigue governs instead, and the seal section stops at 45 mm. On a real design, that number
comes from the bearing maker.

### 5.5 Critical speed
A critical speed is where the rotor's natural frequency equals the running speed, and vibration
builds. This rotor is "stiff": the first critical is at 21,100 rpm, **5.9×** running speed (the rule
used is ≥ 1.2×). These are **dry** criticals: in service, the wear rings act like extra fluid bearings
(the Lomakin effect), which only raises them, so dry is the conservative case. The impeller mass is an
estimate (3.5 kg); a heavier casting would lower the critical, but not by anything near 5×.

### 5.6 Fatigue (DE-Goodman, from Shigley)
**Why the bending is fully reversed:** the radial load points in a fixed direction (toward the low
pressure side of the volute), but the shaft turns. Any point on the surface goes from tension to
compression once per revolution. So bending is purely alternating, and torque is steady.

DE-Goodman combines them: 1/n = 16/(πd³)·[2·Kf·Ma/Se + √3·Kfs·Tm/Sut].
- **Se**, the endurance limit, starts at 0.5·Sut and is cut down by **Marin factors**: machined
  surface (rougher means earlier cracks), size (bigger parts are more likely to contain a flaw), and 99%
  reliability (0.814).
- **Kf, Kfs** are stress-concentration factors: keyseat 2.14/3.0, sharp bearing shoulder 2.7/2.2,
  3 mm fillet 1.7/1.5. Notch sensitivity is taken as 1, which is conservative.
- **Why the bearing shoulder governs:** the bending moment peaks at the inboard bearing (the load times
  the whole overhang), and a bearing shoulder needs a small fillet to clear the bearing's corner, so
  its stress concentration is high. Result: n = 2.4 there and 2.7 at the keyseat.

### 5.7 Bearings
L10 life is the life 90% of bearings reach: L10 = (C/P)³ × 10⁶ revolutions, for ball bearings.
Solving for C gives the catalog rating to look for. API 610 asks for 25,000 h at rated conditions
and 16,000 h at maximum loads. The outboard pair takes the thrust. For back-to-back 40° angular-contact
pairs, the equivalent load is P = Fr + 0.55·Fa (for Fa/Fr ≤ 1.14) or 0.57·Fr + 0.93·Fa. The tool reports
the C required (about 49 kN each position) rather than inventing a catalog part.

### 5.8 Key
Square key per DIN 6885 for a 30 mm shaft: 8 × 7 mm. Two ways it fails: shear across its width
(τ = 2T/(d·b·L)), and crushing on the hub side (σ = 2T/(d·(h − t1)·L)). Safety factor 4.8.

---

## 6. 3D impeller stress (`pump/impeller_fea.py`)

### 6.1 Geometry
Built in gmsh (OpenCASCADE) from the hydraulic design: flat back shroud, conical front shroud,
wear-ring lip, hub, and 5 blades from the same camber lines as the plan view. The camber line comes
from integrating dθ = dr/(r·tan β), with β varying linearly from inlet to outlet. Simplifications:
cylindrical blades (the same profile at every height), constant thickness, a conical front shroud,
no fillets.

### 6.2 Loads, and why each one
1. **Centrifugal** (ρω²r): the dominant load on shrouds.
2. **Passage pressure** on every wetted surface. It is from rothalpy, the quantity conserved along a
   streamline in a rotating frame: p − p_s = ρ/2·[(u² − u1²) − (w² − w1²)]. Pressure rises
   outward because of the centrifugal term (u²), and because the relative flow slows down (w falls).
3. **Side-gap pressure** on the shroud backs: the outlet pressure minus the drop from fluid spinning
   at half speed. Near the eye this pushes the shrouds inward, because the outside pressure is higher
   than the passage pressure there.
4. **Blade loading**: a blade can only turn the flow if its two faces see different pressures. The
   circulation each blade carries (2π·r2·cu2/z), spread along its length, gives the mean velocity
   difference Δw across it, and the pressure difference ρ·w·Δw ≈ 0.4 MPa.

### 6.3 Solver choices
- **Quadratic (P2) tetrahedra:** straight-edge linear tetrahedra are far too stiff in bending, and
  blades and shrouds are thin bending members.
- **Clamped bore:** simple and slightly over-stiff locally. Hub stresses near the bore are reported
  separately.
- **Iterative solver with multigrid, given the rigid-body modes:** a direct solver for 783k unknowns
  needs a lot of memory. Multigrid needs to know which motions cost no strain energy (3 translations,
  3 rotations) to coarsen correctly. Giving it those cut the solve time 7× (255 s → 37 s on the test
  disc).

### 6.4 Verification and convergence, and why the peaks don't count
- **Rotating disc, clamped bore:** hoop stress within 0.3% of the closed-form solution at every radius.
- **Mesh convergence:** coarse → fine changed the 99th-percentile stresses by 4–5%: **converged**.
  The corner peaks rose 213 → 273 MPa: **not converged**. At a sharp re-entrant corner the
  elasticity solution is mathematically infinite, so the finer the mesh, the higher the number. That
  is called a singularity. Real castings have fillets there. So: use the converged values for the
  margin (79 MPa vs 550 MPa yield, about 7×), and read the corner peaks as "this corner needs a fillet".
- **Mass check:** the solid weighs 3.57 kg, the hand estimate said 3.49 kg.
- **Wear-ring lip growth:** 2–4 µm against 200 µm radial clearance. Spin and pressure will not close
  the clearance.

---

## 7. Where I chose not to model something, and why

| Not modeled | Why | Effect of leaving it out |
|---|---|---|
| CFD | hours per design on this machine; can't sweep 1,000 designs | absolute numbers within about 10–20%; trends more reliable |
| Low-flow recirculation | no 1D model for it | below 10% of rated flow the curve isn't judged |
| Pressure pulsation / impeller fatigue | needs the pulsation amplitude, which 1D doesn't give | the impeller's alternating stress is unknown |
| Wet criticals (Lomakin) | needs ring-seal dynamics | dry criticals are conservative |
| Fillets in 3D | simplicity | corner peaks are singular, see 6.4 |
| Gyroscopic effects | small for a light, stiff rotor | criticals slightly off |

---

## 8. Problems found along the way (each one is an interview story)

1. **Efficiency came out at 86%, too high.** I didn't believe it, since a pump this size is about 80%.
   That led to calibrating against published statistics instead of trusting textbook coefficients.
2. **BEP drifted to 117% of rated flow** after calibration. Real designers move the BEP with the
   volute throat, so the sizing now chooses the throat to put rated flow at BEP.
3. **The inlet was running at −3° incidence when +3° was asked.** The blade angle was set from a flow
   angle that ignored blockage and leakage. Fixed with an iteration on the real flow angle, and a test
   now guards it.
4. **The optimizer chose 4 blades.** Investigated, tried a blade-loading loss, found it made things
   worse, bounded blade count by practice and documented why (section 4).
5. **The stability check failed designs on 0.006% humps**, which is numerical noise. Added a 0.05%
   tolerance and documented it, rather than loosening the check silently.
6. **The hub was a placeholder (31 mm)** when the real shaft needs a 48 mm hub. Fixing it cut 14% of
   the eye area, dropped the baseline's NPSH margin below the limit, and showed the eye-sizing rule
   ignored the hub. The rule was corrected, and the optimizer's gain fell from 1.7 to 1.0 points.
   **Lesson: an early placeholder can quietly decide the answer.**
7. **The 3D mesher failed**: overlapping faces where blade tips met the shroud rim. Fixed by building
   everything oversize and turning the OD in one cut, like a lathe.
8. **The 3D solve was slow** (255 s on a test disc). Giving multigrid the rigid-body modes made it 37 s.

---

## 9. Self-quiz

Cover the answers and say them out loud.

1. **Why does the impeller's head fall as flow rises?** More flow raises cm2. With backswept blades
   that lowers cu2, so Euler head u2·cu2/g falls. Losses that grow with flow add to it.
2. **What is slip and what does it do here?** The flow lags the blade at the outlet because of a
   relative eddy. σ = 0.78 here, so about 22% of the ideal swirl is lost. More blades reduce it.
3. **Why is there a best efficiency point?** Friction grows with Q², shock losses grow away from the
   design flow, and disk friction is a fixed tax. The peak is where they balance.
4. **Why does a bigger eye not always lower NPSHr?** It lowers the inlet velocity but raises the blade
   speed at the eye, so the leading-edge velocity w1 grows. There's an optimum.
5. **What did the optimizer trade, and why did it stop where it did?** It shrank the eye for less leakage
   and friction, which raised NPSHr until the 1 m margin was used up.
6. **Why calibrate only friction losses?** The data is at BEP, where shock losses are near zero, so it
   says nothing about them.
7. **How do you know the calibration is right?** The overall-efficiency correlation, which wasn't used to tune,
   agrees: 81.8% vs 82.0%.
8. **Why does a steep outlet angle cause a hump?** It flattens the ideal head line, and the low-flow
   losses then bend the curve over. A humped pump can hunt between two flows.
9. **Why is the radial load near zero at BEP?** The volute is designed for that flow, so pressure is
   even around the impeller. Off-design, it is uneven.
10. **Why size the shaft with the maximum impeller?** API 610 requires the pump to take a 5%-more-head
    impeller later, and the shaft must survive it.
11. **Why is shaft bending fully reversed?** The load's direction is fixed and the shaft rotates, so each
    point alternates tension and compression every revolution.
12. **Where does the shaft fail first in fatigue, and why there?** The inboard bearing shoulder:
    the bending moment peaks there, and a bearing shoulder needs a small fillet, so stress concentration is high.
13. **What sets the seal-section diameter?** Deflection at the seal faces (50 µm) for soft bearings;
    fatigue at the bearing shoulder once the bearings are stiff (3×10⁸ N/m and up).
14. **Why is the thrust bearing sized for 2.6 kN when nominal thrust is 0.16 kN?** The residual thrust on a
    balanced impeller is uncertain, so it's sized for the ±10% shroud-force allowance.
15. **Dry vs wet critical speed?** Dry has no fluid in the wear rings. In service, the rings add stiffness
    (Lomakin), which raises the critical, so dry is conservative.
16. **Why quadratic elements in the 3D model?** Linear tetrahedra are far too stiff in bending, and the
    blades and shrouds bend.
17. **The peak stress was 313 MPa. Is that real?** No. It rose with mesh refinement, which marks a
    singularity at a sharp corner. The converged value is 79 MPa, and the corner needs a fillet.
18. **How did you verify the 3D solver?** A rotating disc with a clamped bore, against the closed-form solution:
    within 0.3%.
19. **What would you do next?** Validate against a published pump test curve; model fillets and
    blade-pass pressure pulsation for impeller fatigue; get real bearing stiffness.
20. **Biggest weakness of the whole study?** It's a 1D hydraulic model with correlations (10–20%
    absolute error), and the drive-end layout is assumed. The trends and the method are the value,
    not the exact numbers.
21. **How did you build it?** I directed an AI coding assistant, then checked the engineering, the tests
    and the results myself. (Say this plainly if asked.)
