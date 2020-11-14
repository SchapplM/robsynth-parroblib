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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2020-08-06 17:43
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3PRRRR8V2G1P3A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(8,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G1P3A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G1P3A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G1P3A1_Jinv: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G1P3A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G1P3A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:43:31
% EndTime: 2020-08-06 17:43:32
% DurationCPUTime: 0.27s
% Computational Cost: add. (144->61), mult. (306->139), div. (6->3), fcn. (318->22), ass. (0->70)
t40 = sin(qJ(2,1));
t46 = cos(qJ(2,1));
t47 = pkin(7) + pkin(6);
t18 = pkin(2) * t46 + t40 * t47;
t28 = sin(pkin(8));
t30 = cos(pkin(8));
t15 = pkin(2) * t40 - t47 * t46;
t31 = cos(pkin(4));
t29 = sin(pkin(4));
t39 = sin(qJ(3,1));
t60 = t29 * t39;
t48 = pkin(3) * t60 - t15 * t31;
t72 = t18 * t30 + t48 * t28;
t38 = sin(qJ(2,2));
t44 = cos(qJ(2,2));
t17 = pkin(2) * t44 + t38 * t47;
t14 = pkin(2) * t38 - t47 * t44;
t37 = sin(qJ(3,2));
t61 = t29 * t37;
t49 = pkin(3) * t61 - t14 * t31;
t71 = t17 * t30 + t49 * t28;
t36 = sin(qJ(2,3));
t42 = cos(qJ(2,3));
t16 = pkin(2) * t42 + t36 * t47;
t13 = pkin(2) * t36 - t47 * t42;
t35 = sin(qJ(3,3));
t62 = t29 * t35;
t50 = pkin(3) * t62 - t13 * t31;
t70 = t16 * t30 + t50 * t28;
t41 = cos(qJ(3,3));
t25 = t41 ^ 2;
t69 = pkin(3) * t25;
t43 = cos(qJ(3,2));
t26 = t43 ^ 2;
t68 = pkin(3) * t26;
t45 = cos(qJ(3,1));
t27 = t45 ^ 2;
t67 = pkin(3) * t27;
t66 = pkin(3) * t29;
t59 = t31 * t35;
t58 = t31 * t36;
t57 = t31 * t37;
t56 = t31 * t38;
t55 = t31 * t39;
t54 = t31 * t40;
t53 = pkin(2) * t62;
t52 = pkin(2) * t61;
t51 = pkin(2) * t60;
t34 = legFrame(1,3);
t33 = legFrame(2,3);
t32 = legFrame(3,3);
t24 = cos(t34);
t23 = cos(t33);
t22 = cos(t32);
t21 = sin(t34);
t20 = sin(t33);
t19 = sin(t32);
t12 = t28 * t46 + t30 * t54;
t11 = t28 * t44 + t30 * t56;
t10 = t28 * t42 + t30 * t58;
t9 = t28 * t54 - t30 * t46;
t8 = t28 * t56 - t30 * t44;
t7 = t28 * t58 - t30 * t42;
t6 = t28 * t18 - t48 * t30;
t5 = t28 * t17 - t49 * t30;
t4 = t28 * t16 - t50 * t30;
t3 = 0.1e1 / (t40 * t27 * t66 + (pkin(3) * t55 + t15 * t29) * t45 + pkin(2) * t55);
t2 = 0.1e1 / (t38 * t26 * t66 + (pkin(3) * t57 + t14 * t29) * t43 + pkin(2) * t57);
t1 = 0.1e1 / (t36 * t25 * t66 + (pkin(3) * t59 + t13 * t29) * t41 + pkin(2) * t59);
t63 = [(-(t21 * t12 + t9 * t24) * t67 + (-t6 * t21 + t24 * t72) * t45 + (t30 * t21 + t24 * t28) * t51) * t3, ((t12 * t24 - t21 * t9) * t67 + (t21 * t72 + t6 * t24) * t45 - (-t28 * t21 + t24 * t30) * t51) * t3, 1; (-(t20 * t11 + t8 * t23) * t68 + (-t5 * t20 + t23 * t71) * t43 + (t30 * t20 + t23 * t28) * t52) * t2, ((t11 * t23 - t20 * t8) * t68 + (t20 * t71 + t5 * t23) * t43 - (-t28 * t20 + t23 * t30) * t52) * t2, 1; (-(t19 * t10 + t7 * t22) * t69 + (-t4 * t19 + t22 * t70) * t41 + (t30 * t19 + t22 * t28) * t53) * t1, ((t10 * t22 - t19 * t7) * t69 + (t19 * t70 + t4 * t22) * t41 - (-t28 * t19 + t22 * t30) * t53) * t1, 1;];
Jinv  = t63;
