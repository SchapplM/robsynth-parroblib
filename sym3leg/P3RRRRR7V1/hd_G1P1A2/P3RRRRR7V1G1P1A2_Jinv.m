% Analytische Jacobi-Matrix für parallelen Roboter
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorpose und aktiven Gelenkkoordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d3,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% Jinv [3x3]
%   Analytische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 03:40
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRRRR7V1G1P1A2_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(5,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V1G1P1A2_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V1G1P1A2_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRRRR7V1G1P1A2_Jinv: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V1G1P1A2_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V1G1P1A2_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:40:12
% EndTime: 2020-08-07 03:40:12
% DurationCPUTime: 0.16s
% Computational Cost: add. (237->97), mult. (156->54), div. (18->7), fcn. (84->69), ass. (0->46)
t52 = -2 * pkin(4) - 2 * pkin(5);
t47 = 0.1e1 / pkin(1);
t51 = t47 / 0.2e1;
t28 = legFrame(3,3) + qJ(1,3);
t29 = legFrame(2,3) + qJ(1,2);
t30 = legFrame(1,3) + qJ(1,1);
t31 = sin(qJ(2,3) + qJ(3,3));
t38 = 0.2e1 * qJ(3,3);
t50 = 0.1e1 / ((sin(qJ(2,3)) - sin(t38 + qJ(2,3))) * pkin(2) + (sin(qJ(2,3) - qJ(3,3)) - t31) * pkin(1)) * t51;
t32 = sin(qJ(2,2) + qJ(3,2));
t41 = 0.2e1 * qJ(3,2);
t49 = 0.1e1 / ((-sin(t41 + qJ(2,2)) + sin(qJ(2,2))) * pkin(2) + (sin(qJ(2,2) - qJ(3,2)) - t32) * pkin(1)) * t51;
t33 = sin(qJ(2,1) + qJ(3,1));
t44 = 0.2e1 * qJ(3,1);
t48 = 0.1e1 / ((sin(qJ(2,1)) - sin(t44 + qJ(2,1))) * pkin(2) + (sin(qJ(2,1) - qJ(3,1)) - t33) * pkin(1)) * t51;
t23 = -qJ(3,3) + t28;
t22 = qJ(3,3) + t28;
t25 = -qJ(3,2) + t29;
t24 = qJ(3,2) + t29;
t27 = -qJ(3,1) + t30;
t26 = qJ(3,1) + t30;
t46 = -0.2e1 * qJ(2,1);
t45 = 0.2e1 * qJ(2,1);
t43 = -0.2e1 * qJ(2,2);
t42 = 0.2e1 * qJ(2,2);
t40 = -0.2e1 * qJ(2,3);
t39 = 0.2e1 * qJ(2,3);
t21 = -0.2e1 * qJ(3,1) + t46 + t30;
t20 = -qJ(2,1) + t27;
t19 = qJ(2,1) + t26;
t18 = t44 + t45 + t30;
t17 = -0.2e1 * qJ(3,2) + t43 + t29;
t16 = -qJ(2,2) + t25;
t15 = qJ(2,2) + t24;
t14 = t41 + t42 + t29;
t13 = -0.2e1 * qJ(3,3) + t40 + t28;
t12 = -qJ(2,3) + t23;
t11 = qJ(2,3) + t22;
t10 = t38 + t39 + t28;
t9 = t46 + t27;
t8 = t45 + t26;
t7 = t43 + t25;
t6 = t42 + t24;
t5 = t40 + t23;
t4 = t39 + t22;
t1 = [((sin(t20) + sin(t19)) * t52 + (-cos(t21) - cos(t18) - 0.2e1 * cos(t30)) * pkin(2) + (-cos(t9) - cos(t8) - cos(t27) - cos(t26)) * pkin(1)) * t48, ((-cos(t20) - cos(t19)) * t52 + (-sin(t21) - sin(t18) - 0.2e1 * sin(t30)) * pkin(2) + (-sin(t9) - sin(t8) - sin(t27) - sin(t26)) * pkin(1)) * t48, t47 * t33 / sin(qJ(3,1)); ((sin(t16) + sin(t15)) * t52 + (-cos(t17) - cos(t14) - 0.2e1 * cos(t29)) * pkin(2) + (-cos(t7) - cos(t6) - cos(t25) - cos(t24)) * pkin(1)) * t49, ((-cos(t16) - cos(t15)) * t52 + (-sin(t17) - sin(t14) - 0.2e1 * sin(t29)) * pkin(2) + (-sin(t7) - sin(t6) - sin(t25) - sin(t24)) * pkin(1)) * t49, t47 * t32 / sin(qJ(3,2)); ((sin(t12) + sin(t11)) * t52 + (-cos(t13) - cos(t10) - 0.2e1 * cos(t28)) * pkin(2) + (-cos(t5) - cos(t4) - cos(t23) - cos(t22)) * pkin(1)) * t50, ((-cos(t12) - cos(t11)) * t52 + (-sin(t13) - sin(t10) - 0.2e1 * sin(t28)) * pkin(2) + (-sin(t5) - sin(t4) - sin(t23) - sin(t22)) * pkin(1)) * t50, t47 * t31 / sin(qJ(3,3));];
Jinv  = t1;
