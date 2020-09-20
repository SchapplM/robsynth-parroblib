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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-08-06 18:51
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RPRRR9V1G1P1A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(7,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G1P1A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G1P1A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G1P1A1_Jinv: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G1P1A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G1P1A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:50:59
% EndTime: 2020-08-06 18:50:59
% DurationCPUTime: 0.11s
% Computational Cost: add. (111->50), mult. (99->72), div. (12->6), fcn. (93->28), ass. (0->41)
t40 = 2 * pkin(1);
t23 = pkin(7) + qJ(3,3);
t8 = cos(t23);
t39 = pkin(3) * t8;
t24 = pkin(7) + qJ(3,2);
t9 = cos(t24);
t38 = pkin(3) * t9;
t37 = -pkin(5) - pkin(6);
t25 = pkin(7) + qJ(3,1);
t10 = cos(t25);
t36 = pkin(3) * t10;
t35 = 2 * pkin(7);
t34 = cos(qJ(1,1));
t33 = cos(qJ(1,2));
t32 = cos(qJ(1,3));
t31 = sin(qJ(1,1));
t30 = sin(qJ(1,2));
t29 = sin(qJ(1,3));
t28 = legFrame(1,3);
t27 = legFrame(2,3);
t26 = legFrame(3,3);
t22 = -qJ(2,1) + t37;
t21 = -qJ(2,2) + t37;
t20 = -qJ(2,3) + t37;
t19 = 0.1e1 / t22;
t18 = 0.1e1 / t21;
t17 = 0.1e1 / t20;
t16 = cos(t28);
t15 = cos(t27);
t14 = cos(t26);
t13 = sin(t28);
t12 = sin(t27);
t11 = sin(t26);
t7 = cos(pkin(7)) * pkin(2) + pkin(1);
t6 = -t31 * t22 + t7 * t34;
t5 = -t30 * t21 + t7 * t33;
t4 = -t29 * t20 + t7 * t32;
t3 = t34 * t22 + t31 * t7;
t2 = t33 * t21 + t30 * t7;
t1 = t32 * t20 + t29 * t7;
t41 = [-((-t13 * t31 + t16 * t34) * t36 + t6 * t16 - t13 * t3) * t19, -((t13 * t34 + t16 * t31) * t36 + t3 * t16 + t13 * t6) * t19, -t19 * (pkin(3) * sin(0.2e1 * t25) + sin(t25) * t40 + (sin((t35 + qJ(3,1))) + sin(qJ(3,1))) * pkin(2)) / t10 / 0.2e1; -((-t12 * t30 + t15 * t33) * t38 + t5 * t15 - t12 * t2) * t18, -((t12 * t33 + t15 * t30) * t38 + t2 * t15 + t12 * t5) * t18, -t18 * (sin(t24) * t40 + pkin(3) * sin(0.2e1 * t24) + (sin((t35 + qJ(3,2))) + sin(qJ(3,2))) * pkin(2)) / t9 / 0.2e1; -((-t11 * t29 + t14 * t32) * t39 + t4 * t14 - t11 * t1) * t17, -((t11 * t32 + t14 * t29) * t39 + t1 * t14 + t11 * t4) * t17, -t17 * (pkin(3) * sin(0.2e1 * t23) + sin(t23) * t40 + (sin((t35 + qJ(3,3))) + sin(qJ(3,3))) * pkin(2)) / t8 / 0.2e1;];
Jinv  = t41;
