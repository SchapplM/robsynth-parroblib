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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2020-08-06 21:11
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRPRR8V2G1P1A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(7,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G1P1A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G1P1A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G1P1A1_Jinv: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G1P1A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G1P1A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:11:13
% EndTime: 2020-08-06 21:11:13
% DurationCPUTime: 0.18s
% Computational Cost: add. (117->58), mult. (132->88), div. (12->6), fcn. (99->33), ass. (0->48)
t50 = -pkin(5) - pkin(6);
t28 = qJ(2,3) + pkin(7);
t49 = pkin(3) * cos(t28);
t29 = qJ(2,2) + pkin(7);
t48 = pkin(3) * cos(t29);
t30 = qJ(2,1) + pkin(7);
t47 = pkin(3) * cos(t30);
t46 = 0.2e1 * pkin(2) * pkin(3);
t45 = 2 * pkin(1);
t44 = pkin(2) ^ 2;
t43 = pkin(3) ^ 2;
t42 = 2 * qJ(2,1);
t41 = 2 * qJ(2,2);
t40 = 2 * qJ(2,3);
t39 = cos(qJ(1,1));
t38 = cos(qJ(1,2));
t37 = cos(qJ(1,3));
t36 = sin(qJ(1,1));
t35 = sin(qJ(1,2));
t34 = sin(qJ(1,3));
t33 = legFrame(1,3);
t32 = legFrame(2,3);
t31 = legFrame(3,3);
t27 = -qJ(3,1) + t50;
t26 = -qJ(3,2) + t50;
t25 = -qJ(3,3) + t50;
t24 = cos(qJ(2,1)) * pkin(2);
t23 = cos(qJ(2,2)) * pkin(2);
t22 = cos(qJ(2,3)) * pkin(2);
t21 = 0.1e1 / t27;
t20 = 0.1e1 / t26;
t19 = 0.1e1 / t25;
t18 = cos(t33);
t17 = cos(t32);
t16 = cos(t31);
t15 = sin(t33);
t14 = sin(t32);
t13 = sin(t31);
t9 = t24 + pkin(1);
t8 = t23 + pkin(1);
t7 = t22 + pkin(1);
t6 = -t36 * t27 + t9 * t39;
t5 = -t35 * t26 + t8 * t38;
t4 = -t34 * t25 + t7 * t37;
t3 = t39 * t27 + t36 * t9;
t2 = t38 * t26 + t35 * t8;
t1 = t37 * t25 + t34 * t7;
t10 = [-((-t15 * t36 + t18 * t39) * t47 + t6 * t18 - t3 * t15) * t21, -((t15 * t39 + t18 * t36) * t47 + t3 * t18 + t6 * t15) * t21, -(sin((t42 + pkin(7))) * t46 + t44 * sin(t42) + t43 * sin(0.2e1 * t30) + (sin(t30) * pkin(3) + sin(qJ(2,1)) * pkin(2)) * t45) * t21 / (t24 + t47) / 0.2e1; -((-t14 * t35 + t17 * t38) * t48 + t5 * t17 - t2 * t14) * t20, -((t14 * t38 + t17 * t35) * t48 + t2 * t17 + t5 * t14) * t20, -(sin((t41 + pkin(7))) * t46 + t44 * sin(t41) + t43 * sin(0.2e1 * t29) + (sin(qJ(2,2)) * pkin(2) + sin(t29) * pkin(3)) * t45) * t20 / (t23 + t48) / 0.2e1; -((-t13 * t34 + t16 * t37) * t49 + t4 * t16 - t1 * t13) * t19, -((t13 * t37 + t16 * t34) * t49 + t1 * t16 + t4 * t13) * t19, -(sin((t40 + pkin(7))) * t46 + t44 * sin(t40) + t43 * sin(0.2e1 * t28) + (sin(qJ(2,3)) * pkin(2) + sin(t28) * pkin(3)) * t45) * t19 / (t22 + t49) / 0.2e1;];
Jinv  = t10;
