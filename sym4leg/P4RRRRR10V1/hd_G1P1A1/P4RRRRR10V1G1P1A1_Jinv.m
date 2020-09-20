% Analytische Jacobi-Matrix für parallelen Roboter
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorpose und aktiven Gelenkkoordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Input:
% xP [4x1]
%   Generalized platform coordinates
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d1,d2,d4]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% Jinv [4x4]
%   Analytische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 13:25
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P4RRRRR10V1G1P1A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(6,1),zeros(4,3),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4RRRRR10V1G1P1A1_Jinv: qJ has to be [3x4] (double)');
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4RRRRR10V1G1P1A1_Jinv: xP has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P4RRRRR10V1G1P1A1_Jinv: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4RRRRR10V1G1P1A1_Jinv: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4RRRRR10V1G1P1A1_Jinv: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 13:23:46
% EndTime: 2020-08-07 13:23:46
% DurationCPUTime: 0.38s
% Computational Cost: add. (228->76), mult. (644->157), div. (20->4), fcn. (620->36), ass. (0->91)
t43 = sin(pkin(3));
t61 = sin(qJ(3,1));
t70 = cos(qJ(3,1));
t71 = cos(qJ(2,1));
t86 = t70 * t71;
t62 = sin(qJ(2,1));
t28 = pkin(2) * t62 * t70 - t71 * pkin(6);
t44 = cos(pkin(3));
t94 = t28 * t44;
t98 = pkin(6) * t62;
t12 = 0.1e1 / (pkin(1) * t94 + (pkin(5) * t98 + (-pkin(1) * t61 + pkin(5) * t86) * pkin(2)) * t43);
t58 = sin(qJ(3,2));
t67 = cos(qJ(3,2));
t68 = cos(qJ(2,2));
t87 = t67 * t68;
t59 = sin(qJ(2,2));
t27 = pkin(2) * t59 * t67 - t68 * pkin(6);
t95 = t27 * t44;
t99 = pkin(6) * t59;
t11 = 0.1e1 / (pkin(1) * t95 + (pkin(5) * t99 + (-pkin(1) * t58 + pkin(5) * t87) * pkin(2)) * t43);
t56 = sin(qJ(2,3));
t100 = pkin(6) * t56;
t55 = sin(qJ(3,3));
t64 = cos(qJ(3,3));
t65 = cos(qJ(2,3));
t88 = t64 * t65;
t26 = pkin(2) * t56 * t64 - t65 * pkin(6);
t96 = t26 * t44;
t10 = 0.1e1 / (pkin(1) * t96 + (pkin(5) * t100 + (-pkin(1) * t55 + pkin(5) * t88) * pkin(2)) * t43);
t50 = sin(qJ(2,4));
t101 = pkin(6) * t50;
t49 = sin(qJ(3,4));
t52 = cos(qJ(3,4));
t53 = cos(qJ(2,4));
t89 = t52 * t53;
t25 = pkin(2) * t50 * t52 - t53 * pkin(6);
t97 = t25 * t44;
t9 = 0.1e1 / (pkin(1) * t97 + (pkin(5) * t101 + (-pkin(1) * t49 + pkin(5) * t89) * pkin(2)) * t43);
t102 = pkin(2) * t44;
t93 = t43 * t49;
t92 = t43 * t55;
t91 = t43 * t58;
t90 = t43 * t61;
t81 = koppelP(1,1);
t80 = koppelP(2,1);
t79 = koppelP(3,1);
t78 = koppelP(4,1);
t77 = koppelP(1,2);
t76 = koppelP(2,2);
t75 = koppelP(3,2);
t74 = koppelP(4,2);
t73 = xP(4);
t72 = cos(qJ(1,1));
t69 = cos(qJ(1,2));
t66 = cos(qJ(1,3));
t63 = sin(qJ(1,1));
t60 = sin(qJ(1,2));
t57 = sin(qJ(1,3));
t54 = cos(qJ(1,4));
t51 = sin(qJ(1,4));
t48 = legFrame(1,3);
t47 = legFrame(2,3);
t46 = legFrame(3,3);
t45 = legFrame(4,3);
t42 = cos(t73);
t41 = sin(t73);
t36 = cos(t48);
t35 = cos(t47);
t34 = cos(t46);
t33 = cos(t45);
t32 = sin(t48);
t31 = sin(t47);
t30 = sin(t46);
t29 = sin(t45);
t24 = t32 * t72 + t36 * t63;
t23 = t31 * t69 + t35 * t60;
t22 = t30 * t66 + t34 * t57;
t21 = t32 * t63 - t36 * t72;
t20 = t31 * t60 - t35 * t69;
t19 = t30 * t57 - t34 * t66;
t18 = t29 * t54 + t33 * t51;
t17 = t29 * t51 - t33 * t54;
t8 = t21 * t98 + t24 * t94 + (t21 * t86 - t24 * t90) * pkin(2);
t7 = t20 * t99 + t23 * t95 + (t20 * t87 - t23 * t91) * pkin(2);
t6 = t19 * t100 + t22 * t96 + (t19 * t88 - t22 * t92) * pkin(2);
t5 = -t24 * t98 + t21 * t94 + (-t21 * t90 - t24 * t86) * pkin(2);
t4 = -t23 * t99 + t20 * t95 + (-t20 * t91 - t23 * t87) * pkin(2);
t3 = -t22 * t100 + t19 * t96 + (-t19 * t92 - t22 * t88) * pkin(2);
t2 = t17 * t101 + t18 * t97 + (t17 * t89 - t18 * t93) * pkin(2);
t1 = -t18 * t101 + t17 * t97 + (-t17 * t93 - t18 * t89) * pkin(2);
t13 = [-t8 * t12, -t5 * t12, (t102 * t61 + t28 * t43) * t12, -(t8 * (-t41 * t81 - t42 * t77) + t5 * (-t41 * t77 + t42 * t81)) * t12; -t7 * t11, -t4 * t11, (t102 * t58 + t27 * t43) * t11, -(t7 * (-t41 * t80 - t42 * t76) + t4 * (-t41 * t76 + t42 * t80)) * t11; -t6 * t10, -t3 * t10, (t102 * t55 + t26 * t43) * t10, -(t6 * (-t41 * t79 - t42 * t75) + t3 * (-t41 * t75 + t42 * t79)) * t10; -t2 * t9, -t1 * t9, (t102 * t49 + t25 * t43) * t9, -(t2 * (-t41 * t78 - t42 * t74) + t1 * (-t41 * t74 + t42 * t78)) * t9;];
Jinv  = t13;
