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
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
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
% Datum: 2020-08-07 00:53
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRRRR10V2G1P1A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(8,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR10V2G1P1A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR10V2G1P1A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3RRRRR10V2G1P1A1_Jinv: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR10V2G1P1A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR10V2G1P1A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 00:51:56
% EndTime: 2020-08-07 00:51:56
% DurationCPUTime: 0.36s
% Computational Cost: add. (234->88), mult. (483->185), div. (9->3), fcn. (465->26), ass. (0->78)
t23 = cos(pkin(4));
t77 = pkin(1) * t23;
t76 = pkin(2) * t23;
t36 = cos(qJ(3,3));
t19 = t36 ^ 2;
t75 = pkin(3) * t19;
t39 = cos(qJ(3,2));
t20 = t39 ^ 2;
t74 = pkin(3) * t20;
t42 = cos(qJ(3,1));
t21 = t42 ^ 2;
t73 = pkin(3) * t21;
t22 = sin(pkin(4));
t72 = pkin(3) * t22;
t71 = pkin(6) * t22;
t27 = sin(qJ(3,3));
t70 = t27 * pkin(3);
t30 = sin(qJ(3,2));
t69 = t30 * pkin(3);
t33 = sin(qJ(3,1));
t68 = t33 * pkin(3);
t37 = cos(qJ(2,3));
t67 = t37 * pkin(2);
t40 = cos(qJ(2,2));
t66 = t40 * pkin(2);
t43 = cos(qJ(2,1));
t65 = t43 * pkin(2);
t64 = t22 * t27;
t63 = t22 * t30;
t62 = t22 * t33;
t61 = t23 * t27;
t28 = sin(qJ(2,3));
t60 = t23 * t28;
t59 = t23 * t30;
t31 = sin(qJ(2,2));
t58 = t23 * t31;
t57 = t23 * t33;
t34 = sin(qJ(2,1));
t56 = t23 * t34;
t45 = pkin(8) + pkin(7);
t55 = t23 * t45;
t54 = t45 * t28;
t53 = t45 * t31;
t52 = t45 * t34;
t24 = legFrame(3,3);
t13 = sin(t24);
t16 = cos(t24);
t29 = sin(qJ(1,3));
t38 = cos(qJ(1,3));
t4 = -t13 * t29 + t16 * t38;
t51 = t4 * t64;
t7 = t38 * t13 + t29 * t16;
t50 = t7 * t64;
t25 = legFrame(2,3);
t14 = sin(t25);
t17 = cos(t25);
t32 = sin(qJ(1,2));
t41 = cos(qJ(1,2));
t5 = -t14 * t32 + t17 * t41;
t49 = t5 * t63;
t8 = t41 * t14 + t32 * t17;
t48 = t8 * t63;
t26 = legFrame(1,3);
t15 = sin(t26);
t18 = cos(t26);
t35 = sin(qJ(1,1));
t44 = cos(qJ(1,1));
t6 = -t15 * t35 + t18 * t44;
t47 = t6 * t62;
t9 = t44 * t15 + t35 * t18;
t46 = t9 * t62;
t12 = -t34 * pkin(2) + t45 * t43;
t11 = -t31 * pkin(2) + t45 * t40;
t10 = -t28 * pkin(2) + t45 * t37;
t3 = 0.1e1 / (-(pkin(1) * t56 + t43 * t71) * t73 + (((-pkin(6) + t68) * t65 - pkin(6) * t52 + pkin(1) * t68) * t22 + t12 * t77) * t42 + (pkin(1) + t52 + t65) * pkin(2) * t62);
t2 = 0.1e1 / (-(pkin(1) * t58 + t40 * t71) * t74 + (((-pkin(6) + t69) * t66 - pkin(6) * t53 + pkin(1) * t69) * t22 + t11 * t77) * t39 + (pkin(1) + t53 + t66) * pkin(2) * t63);
t1 = 0.1e1 / (-(pkin(1) * t60 + t37 * t71) * t75 + (((-pkin(6) + t70) * t67 - pkin(6) * t54 + pkin(1) * t70) * t22 + t10 * t77) * t36 + (pkin(1) + t54 + t67) * pkin(2) * t64);
t78 = [(-(t6 * t43 - t9 * t56) * t73 + (-pkin(3) * t46 + (-pkin(2) * t6 - t9 * t55) * t43 + t34 * (-t45 * t6 + t9 * t76)) * t42 - pkin(2) * t46) * t3, (-(t9 * t43 + t6 * t56) * t73 + (pkin(3) * t47 + (-pkin(2) * t9 + t6 * t55) * t43 - (t45 * t9 + t6 * t76) * t34) * t42 + pkin(2) * t47) * t3, (-t34 * t21 * t72 + (-pkin(3) * t57 + t22 * t12) * t42 - pkin(2) * t57) * t3; (-(t5 * t40 - t8 * t58) * t74 + (-pkin(3) * t48 + (-pkin(2) * t5 - t8 * t55) * t40 + t31 * (-t45 * t5 + t8 * t76)) * t39 - pkin(2) * t48) * t2, (-(t8 * t40 + t5 * t58) * t74 + (pkin(3) * t49 + (-pkin(2) * t8 + t5 * t55) * t40 - (t45 * t8 + t5 * t76) * t31) * t39 + pkin(2) * t49) * t2, (-t31 * t20 * t72 + (-pkin(3) * t59 + t22 * t11) * t39 - pkin(2) * t59) * t2; (-(t4 * t37 - t7 * t60) * t75 + (-pkin(3) * t50 + (-pkin(2) * t4 - t7 * t55) * t37 + t28 * (-t45 * t4 + t7 * t76)) * t36 - pkin(2) * t50) * t1, (-(t7 * t37 + t4 * t60) * t75 + (pkin(3) * t51 + (-pkin(2) * t7 + t4 * t55) * t37 - (t4 * t76 + t45 * t7) * t28) * t36 + pkin(2) * t51) * t1, (-t28 * t19 * t72 + (-pkin(3) * t61 + t22 * t10) * t36 - pkin(2) * t61) * t1;];
Jinv  = t78;
