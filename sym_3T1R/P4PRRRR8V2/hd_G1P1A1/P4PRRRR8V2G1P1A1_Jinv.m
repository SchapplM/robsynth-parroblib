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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2020-08-07 11:15
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P4PRRRR8V2G1P1A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(8,1),zeros(4,3),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G1P1A1_Jinv: qJ has to be [3x4] (double)');
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G1P1A1_Jinv: xP has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G1P1A1_Jinv: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G1P1A1_Jinv: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G1P1A1_Jinv: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:15:25
% EndTime: 2020-08-07 11:15:25
% DurationCPUTime: 0.73s
% Computational Cost: add. (364->137), mult. (760->302), div. (16->8), fcn. (832->30), ass. (0->135)
t72 = sin(qJ(2,1));
t78 = cos(qJ(2,1));
t79 = pkin(7) + pkin(6);
t36 = pkin(2) * t78 + t72 * t79;
t55 = sin(pkin(8));
t57 = cos(pkin(8));
t56 = sin(pkin(4));
t71 = sin(qJ(3,1));
t126 = t56 * t71;
t105 = t79 * t78;
t33 = pkin(2) * t72 - t105;
t58 = cos(pkin(4));
t89 = pkin(3) * t126 - t33 * t58;
t142 = t36 * t57 + t55 * t89;
t70 = sin(qJ(2,2));
t76 = cos(qJ(2,2));
t35 = pkin(2) * t76 + t70 * t79;
t69 = sin(qJ(3,2));
t127 = t56 * t69;
t106 = t79 * t76;
t32 = pkin(2) * t70 - t106;
t90 = pkin(3) * t127 - t32 * t58;
t141 = t35 * t57 + t55 * t90;
t68 = sin(qJ(2,3));
t74 = cos(qJ(2,3));
t34 = pkin(2) * t74 + t68 * t79;
t67 = sin(qJ(3,3));
t128 = t56 * t67;
t107 = t79 * t74;
t31 = pkin(2) * t68 - t107;
t91 = pkin(3) * t128 - t31 * t58;
t140 = t34 * t57 + t55 * t91;
t64 = sin(qJ(2,4));
t66 = cos(qJ(2,4));
t30 = pkin(2) * t66 + t64 * t79;
t63 = sin(qJ(3,4));
t129 = t56 * t63;
t108 = t79 * t66;
t29 = pkin(2) * t64 - t108;
t92 = pkin(3) * t129 - t29 * t58;
t139 = t30 * t57 + t55 * t92;
t65 = cos(qJ(3,4));
t51 = t65 ^ 2;
t138 = pkin(3) * t51;
t73 = cos(qJ(3,3));
t52 = t73 ^ 2;
t137 = pkin(3) * t52;
t75 = cos(qJ(3,2));
t53 = t75 ^ 2;
t136 = pkin(3) * t53;
t77 = cos(qJ(3,1));
t54 = t77 ^ 2;
t135 = pkin(3) * t54;
t134 = pkin(3) * t56;
t125 = t58 * t63;
t124 = t58 * t64;
t123 = t58 * t67;
t122 = t58 * t68;
t121 = t58 * t69;
t120 = t58 * t70;
t119 = t58 * t71;
t118 = t58 * t72;
t117 = t58 * t79;
t116 = t64 * t65;
t115 = t65 * t66;
t114 = t68 * t73;
t113 = t70 * t75;
t112 = t72 * t77;
t111 = t73 * t74;
t110 = t75 * t76;
t109 = t77 * t78;
t104 = pkin(2) * t129;
t103 = pkin(2) * t128;
t102 = pkin(2) * t127;
t101 = pkin(2) * t126;
t37 = t65 * pkin(3) + pkin(2);
t100 = t37 * t129;
t38 = t73 * pkin(3) + pkin(2);
t99 = t38 * t128;
t39 = t75 * pkin(3) + pkin(2);
t98 = t39 * t127;
t40 = t77 * pkin(3) + pkin(2);
t97 = t40 * t126;
t96 = t37 * t58;
t95 = t38 * t58;
t94 = t39 * t58;
t93 = t40 * t58;
t88 = koppelP(1,1);
t87 = koppelP(2,1);
t86 = koppelP(3,1);
t85 = koppelP(4,1);
t84 = koppelP(1,2);
t83 = koppelP(2,2);
t82 = koppelP(3,2);
t81 = koppelP(4,2);
t80 = xP(4);
t62 = legFrame(1,3);
t61 = legFrame(2,3);
t60 = legFrame(3,3);
t59 = legFrame(4,3);
t50 = cos(t80);
t49 = sin(t80);
t48 = cos(t62);
t47 = cos(t61);
t46 = cos(t60);
t45 = cos(t59);
t44 = sin(t62);
t43 = sin(t61);
t42 = sin(t60);
t41 = sin(t59);
t28 = t118 * t57 + t55 * t78;
t27 = t120 * t57 + t55 * t76;
t26 = t122 * t57 + t55 * t74;
t25 = t118 * t55 - t57 * t78;
t24 = t120 * t55 - t57 * t76;
t23 = t122 * t55 - t57 * t74;
t22 = t124 * t57 + t55 * t66;
t21 = t124 * t55 - t57 * t66;
t20 = t57 * t44 + t48 * t55;
t19 = t57 * t43 + t47 * t55;
t18 = t57 * t42 + t46 * t55;
t17 = t57 * t41 + t45 * t55;
t16 = -t55 * t44 + t48 * t57;
t15 = -t55 * t43 + t47 * t57;
t14 = -t55 * t42 + t46 * t57;
t13 = -t55 * t41 + t45 * t57;
t8 = t36 * t55 - t57 * t89;
t7 = t35 * t55 - t57 * t90;
t6 = t34 * t55 - t57 * t91;
t5 = t30 * t55 - t57 * t92;
t4 = 0.1e1 / (t72 * t54 * t134 + (pkin(3) * t119 + t33 * t56) * t77 + pkin(2) * t119);
t3 = 0.1e1 / (t70 * t53 * t134 + (pkin(3) * t121 + t32 * t56) * t75 + pkin(2) * t121);
t2 = 0.1e1 / (t68 * t52 * t134 + (pkin(3) * t123 + t31 * t56) * t73 + pkin(2) * t123);
t1 = 0.1e1 / (t64 * t51 * t134 + (pkin(3) * t125 + t29 * t56) * t65 + pkin(2) * t125);
t9 = [(-(t25 * t48 + t44 * t28) * t135 + (t142 * t48 - t8 * t44) * t77 + t20 * t101) * t4, ((-t44 * t25 + t28 * t48) * t135 + (t142 * t44 + t8 * t48) * t77 - t16 * t101) * t4, 1, (-((t117 * t20 + t16 * t40) * t109 - (-t16 * t79 + t20 * t93) * t112 + t20 * t97) * (t49 * t88 + t50 * t84) + ((-t117 * t16 + t20 * t40) * t109 + (t16 * t93 + t20 * t79) * t112 - t16 * t97) * (-t49 * t84 + t50 * t88)) / (-t77 * t56 * t105 + (t112 * t56 + t119) * t40); (-(t24 * t47 + t43 * t27) * t136 + (t141 * t47 - t7 * t43) * t75 + t19 * t102) * t3, ((-t43 * t24 + t27 * t47) * t136 + (t141 * t43 + t7 * t47) * t75 - t15 * t102) * t3, 1, (-((t19 * t117 + t39 * t15) * t110 - (-t15 * t79 + t19 * t94) * t113 + t19 * t98) * (t49 * t87 + t50 * t83) + ((-t117 * t15 + t19 * t39) * t110 + (t15 * t94 + t19 * t79) * t113 - t15 * t98) * (-t49 * t83 + t50 * t87)) / (-t75 * t56 * t106 + (t113 * t56 + t121) * t39); (-(t23 * t46 + t42 * t26) * t137 + (t140 * t46 - t6 * t42) * t73 + t18 * t103) * t2, ((-t42 * t23 + t26 * t46) * t137 + (t140 * t42 + t6 * t46) * t73 - t14 * t103) * t2, 1, (-((t117 * t18 + t14 * t38) * t111 - (-t14 * t79 + t18 * t95) * t114 + t18 * t99) * (t49 * t86 + t50 * t82) + ((-t117 * t14 + t18 * t38) * t111 + (t14 * t95 + t18 * t79) * t114 - t14 * t99) * (-t49 * t82 + t50 * t86)) / (-t73 * t56 * t107 + (t114 * t56 + t123) * t38); (-(t21 * t45 + t41 * t22) * t138 + (t139 * t45 - t5 * t41) * t65 + t17 * t104) * t1, ((-t41 * t21 + t22 * t45) * t138 + (t139 * t41 + t5 * t45) * t65 - t13 * t104) * t1, 1, (-((t117 * t17 + t13 * t37) * t115 - (-t13 * t79 + t17 * t96) * t116 + t17 * t100) * (t49 * t85 + t50 * t81) + ((-t117 * t13 + t17 * t37) * t115 + (t13 * t96 + t17 * t79) * t116 - t13 * t100) * (-t49 * t81 + t50 * t85)) / (-t65 * t56 * t108 + (t116 * t56 + t125) * t37);];
Jinv  = t9;
