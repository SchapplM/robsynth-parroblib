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
% Datum: 2020-08-06 18:13
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3PRRRR8V2G3P3A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(8,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G3P3A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G3P3A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G3P3A1_Jinv: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G3P3A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G3P3A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:13:38
% EndTime: 2020-08-06 18:13:38
% DurationCPUTime: 0.25s
% Computational Cost: add. (162->67), mult. (360->147), div. (9->3), fcn. (363->22), ass. (0->66)
t35 = sin(qJ(3,3));
t68 = pkin(2) * t35;
t37 = sin(qJ(3,2));
t67 = pkin(2) * t37;
t39 = sin(qJ(3,1));
t66 = pkin(2) * t39;
t41 = cos(qJ(3,3));
t65 = pkin(3) * t41 ^ 2;
t43 = cos(qJ(3,2));
t64 = pkin(3) * t43 ^ 2;
t45 = cos(qJ(3,1));
t63 = pkin(3) * t45 ^ 2;
t29 = sin(pkin(4));
t62 = pkin(3) * t29;
t28 = sin(pkin(8));
t61 = t28 * t29;
t36 = sin(qJ(2,3));
t60 = t29 * t36;
t38 = sin(qJ(2,2));
t59 = t29 * t38;
t40 = sin(qJ(2,1));
t58 = t29 * t40;
t31 = cos(pkin(4));
t57 = t31 * t35;
t56 = t31 * t36;
t55 = t31 * t37;
t54 = t31 * t38;
t53 = t31 * t39;
t52 = t31 * t40;
t30 = cos(pkin(8));
t51 = pkin(2) * t29 * t30;
t42 = cos(qJ(2,3));
t47 = pkin(7) + pkin(6);
t13 = pkin(2) * t36 - t47 * t42;
t50 = -t13 * t31 + t35 * t62;
t44 = cos(qJ(2,2));
t14 = pkin(2) * t38 - t47 * t44;
t49 = -t14 * t31 + t37 * t62;
t46 = cos(qJ(2,1));
t15 = pkin(2) * t40 - t47 * t46;
t48 = -t15 * t31 + t39 * t62;
t34 = legFrame(1,2);
t33 = legFrame(2,2);
t32 = legFrame(3,2);
t24 = cos(t34);
t23 = cos(t33);
t22 = cos(t32);
t21 = sin(t34);
t20 = sin(t33);
t19 = sin(t32);
t18 = pkin(2) * t46 + t40 * t47;
t17 = pkin(2) * t44 + t38 * t47;
t16 = pkin(2) * t42 + t36 * t47;
t12 = t28 * t52 - t30 * t46;
t11 = t28 * t54 - t30 * t44;
t10 = t28 * t56 - t30 * t42;
t9 = pkin(3) * t53 + t29 * t15;
t8 = pkin(3) * t55 + t29 * t14;
t7 = pkin(3) * t57 + t29 * t13;
t6 = t30 * t18 + t48 * t28;
t5 = t30 * t17 + t49 * t28;
t4 = t30 * t16 + t50 * t28;
t3 = 0.1e1 / (pkin(2) * t53 + t9 * t45 + t58 * t63);
t2 = 0.1e1 / (pkin(2) * t55 + t8 * t43 + t59 * t64);
t1 = 0.1e1 / (pkin(2) * t57 + t7 * t41 + t60 * t65);
t25 = [(-(t12 * t24 - t21 * t58) * t63 + (t9 * t21 + t6 * t24) * t45 + (t31 * t21 + t24 * t61) * t66) * t3, ((t12 * t21 + t24 * t58) * t63 + (-t6 * t21 + t24 * t9) * t45 + (-t21 * t61 + t24 * t31) * t66) * t3, (-(t28 * t46 + t30 * t52) * t63 + (-t18 * t28 + t48 * t30) * t45 + t39 * t51) * t3; (-(t11 * t23 - t20 * t59) * t64 + (t8 * t20 + t5 * t23) * t43 + (t31 * t20 + t23 * t61) * t67) * t2, ((t11 * t20 + t23 * t59) * t64 + (-t5 * t20 + t23 * t8) * t43 + (-t20 * t61 + t23 * t31) * t67) * t2, (-(t28 * t44 + t30 * t54) * t64 + (-t17 * t28 + t49 * t30) * t43 + t37 * t51) * t2; (-(t10 * t22 - t19 * t60) * t65 + (t7 * t19 + t4 * t22) * t41 + (t31 * t19 + t22 * t61) * t68) * t1, ((t10 * t19 + t22 * t60) * t65 + (-t4 * t19 + t22 * t7) * t41 + (-t19 * t61 + t22 * t31) * t68) * t1, (-(t28 * t42 + t30 * t56) * t65 + (-t16 * t28 + t50 * t30) * t41 + t35 * t51) * t1;];
Jinv  = t25;
