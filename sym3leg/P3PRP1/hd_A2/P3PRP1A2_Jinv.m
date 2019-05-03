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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
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
% Datum: 2019-05-03 14:46
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3PRP1A2_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRP1A2_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRP1A2_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PRP1A2_Jinv: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRP1A2_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRP1A2_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:46:01
% EndTime: 2019-05-03 14:46:01
% DurationCPUTime: 0.14s
% Computational Cost: add. (99->52), mult. (198->113), div. (9->3), fcn. (99->14), ass. (0->45)
t45 = 2 * pkin(2);
t44 = pkin(2) ^ 2 + 1;
t24 = cos(qJ(2,3));
t43 = t24 * qJ(3,3);
t25 = cos(qJ(2,2));
t42 = t25 * qJ(3,2);
t26 = cos(qJ(2,1));
t41 = t26 * qJ(3,1);
t40 = -0.2e1 * t43;
t39 = -0.2e1 * t42;
t38 = -0.2e1 * t41;
t36 = koppelP(1,1);
t35 = koppelP(2,1);
t34 = koppelP(3,1);
t33 = koppelP(1,2);
t32 = koppelP(2,2);
t31 = koppelP(3,2);
t30 = qJ(3,1) ^ 2;
t29 = qJ(3,2) ^ 2;
t28 = qJ(3,3) ^ 2;
t27 = xP(3);
t23 = sin(qJ(2,1));
t22 = sin(qJ(2,2));
t21 = sin(qJ(2,3));
t20 = legFrame(1,3);
t19 = legFrame(2,3);
t18 = legFrame(3,3);
t17 = cos(t27);
t16 = sin(t27);
t15 = cos(t20);
t14 = cos(t19);
t13 = cos(t18);
t12 = sin(t20);
t11 = sin(t19);
t10 = sin(t18);
t9 = (pkin(2) * t33) + qJ(3,1) * t36;
t8 = pkin(2) * t36 - qJ(3,1) * t33;
t7 = (pkin(2) * t32) + qJ(3,2) * t35;
t6 = pkin(2) * t35 - qJ(3,2) * t32;
t5 = (pkin(2) * t31) + qJ(3,3) * t34;
t4 = pkin(2) * t34 - qJ(3,3) * t31;
t3 = 0.1e1 / (-t30 + (t23 * qJ(3,1) * t45 + (-t30 + t44) * t26) * t26 - t44);
t2 = 0.1e1 / (-t29 + (t22 * qJ(3,2) * t45 + (-t29 + t44) * t25) * t25 - t44);
t1 = 0.1e1 / (-t28 + (t21 * qJ(3,3) * t45 + (-t28 + t44) * t24) * t24 - t44);
t37 = [(t15 * t38 + t23 * (pkin(2) * t15 + qJ(3,1) * t12)) * t3, (t12 * t38 + t23 * (pkin(2) * t12 - qJ(3,1) * t15)) * t3, (((-t16 * t8 - t9 * t17) * t15 + t12 * (-t16 * t9 + t8 * t17)) * t23 + 0.2e1 * ((t16 * t36 + t17 * t33) * t15 - t12 * (-t16 * t33 + t17 * t36)) * t41) * t3; (t14 * t39 + t22 * (pkin(2) * t14 + qJ(3,2) * t11)) * t2, (t11 * t39 + t22 * (pkin(2) * t11 - qJ(3,2) * t14)) * t2, (((-t16 * t6 - t7 * t17) * t14 + t11 * (-t16 * t7 + t6 * t17)) * t22 + 0.2e1 * ((t16 * t35 + t17 * t32) * t14 - t11 * (-t16 * t32 + t17 * t35)) * t42) * t2; (t13 * t40 + t21 * (pkin(2) * t13 + qJ(3,3) * t10)) * t1, (t10 * t40 + t21 * (pkin(2) * t10 - qJ(3,3) * t13)) * t1, (((-t16 * t4 - t5 * t17) * t13 + t10 * (-t16 * t5 + t4 * t17)) * t21 + 0.2e1 * ((t16 * t34 + t17 * t31) * t13 - t10 * (-t16 * t31 + t17 * t34)) * t43) * t1;];
Jinv  = t37;
