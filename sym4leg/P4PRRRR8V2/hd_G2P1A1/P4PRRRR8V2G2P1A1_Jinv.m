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
% Datum: 2020-08-07 11:24
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P4PRRRR8V2G2P1A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(8,1),zeros(4,3),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G2P1A1_Jinv: qJ has to be [3x4] (double)');
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G2P1A1_Jinv: xP has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G2P1A1_Jinv: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G2P1A1_Jinv: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G2P1A1_Jinv: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:24:22
% EndTime: 2020-08-07 11:24:23
% DurationCPUTime: 0.94s
% Computational Cost: add. (388->157), mult. (864->328), div. (20->8), fcn. (854->30), ass. (0->134)
t57 = sin(qJ(2,4));
t59 = cos(qJ(2,4));
t76 = pkin(7) + pkin(6);
t153 = pkin(2) * t57 - t59 * t76;
t65 = sin(qJ(2,3));
t71 = cos(qJ(2,3));
t152 = pkin(2) * t65 - t71 * t76;
t67 = sin(qJ(2,2));
t73 = cos(qJ(2,2));
t151 = pkin(2) * t67 - t73 * t76;
t72 = cos(qJ(3,2));
t140 = pkin(3) * t72 ^ 2;
t109 = t67 * t140;
t150 = t151 * t72 + t109;
t70 = cos(qJ(3,3));
t141 = pkin(3) * t70 ^ 2;
t110 = t65 * t141;
t149 = t152 * t70 + t110;
t58 = cos(qJ(3,4));
t142 = pkin(3) * t58 ^ 2;
t111 = t57 * t142;
t148 = t153 * t58 + t111;
t52 = sin(pkin(8));
t147 = pkin(2) * t52;
t56 = sin(qJ(3,4));
t146 = pkin(2) * t56;
t64 = sin(qJ(3,3));
t145 = pkin(2) * t64;
t66 = sin(qJ(3,2));
t144 = pkin(2) * t66;
t68 = sin(qJ(3,1));
t143 = pkin(2) * t68;
t74 = cos(qJ(3,1));
t139 = pkin(3) * t74 ^ 2;
t53 = sin(pkin(4));
t138 = pkin(3) * t53;
t137 = pkin(3) * t58;
t136 = pkin(3) * t70;
t135 = pkin(3) * t72;
t54 = cos(pkin(8));
t134 = t53 * t54;
t133 = t53 * t57;
t132 = t53 * t65;
t131 = t53 * t67;
t130 = t53 * t68;
t69 = sin(qJ(2,1));
t129 = t53 * t69;
t55 = cos(pkin(4));
t128 = t54 * t55;
t127 = t55 * t57;
t126 = t55 * t65;
t125 = t55 * t67;
t124 = t55 * t68;
t123 = t55 * t69;
t122 = t56 * t55;
t120 = t64 * t55;
t119 = t66 * t55;
t75 = cos(qJ(2,1));
t116 = t75 * t76;
t115 = t53 * t147;
t19 = t54 * t127 + t52 * t59;
t114 = t19 * t142;
t20 = t54 * t126 + t52 * t71;
t113 = t20 * t141;
t21 = t54 * t125 + t52 * t73;
t112 = t21 * t140;
t60 = legFrame(4,2);
t36 = sin(t60);
t108 = t36 * t122;
t61 = legFrame(3,2);
t37 = sin(t61);
t107 = t37 * t120;
t62 = legFrame(2,2);
t38 = sin(t62);
t106 = t38 * t119;
t40 = cos(t60);
t105 = t40 * t122;
t41 = cos(t61);
t104 = t41 * t120;
t42 = cos(t62);
t103 = t42 * t119;
t102 = t76 * t128;
t101 = pkin(2) + t137;
t100 = pkin(2) + t136;
t99 = pkin(2) + t135;
t23 = -t102 + t147;
t33 = t52 * t76;
t24 = pkin(2) * t128 + t33;
t95 = t23 * t59 + t24 * t57;
t94 = t23 * t71 + t24 * t65;
t93 = t23 * t73 + t24 * t67;
t92 = t101 * t56 * t54;
t91 = t100 * t64 * t54;
t90 = t99 * t66 * t54;
t89 = t56 * t138 - t153 * t55;
t88 = t64 * t138 - t152 * t55;
t87 = t66 * t138 - t151 * t55;
t29 = pkin(2) * t69 - t116;
t86 = pkin(3) * t130 - t29 * t55;
t85 = koppelP(1,1);
t84 = koppelP(2,1);
t83 = koppelP(3,1);
t82 = koppelP(4,1);
t81 = koppelP(1,2);
t80 = koppelP(2,2);
t79 = koppelP(3,2);
t78 = koppelP(4,2);
t77 = xP(4);
t63 = legFrame(1,2);
t47 = cos(t77);
t46 = sin(t77);
t43 = cos(t63);
t39 = sin(t63);
t34 = pkin(3) * t74 + pkin(2);
t32 = pkin(2) * t75 + t69 * t76;
t31 = pkin(2) * t73 + t67 * t76;
t30 = pkin(2) * t71 + t65 * t76;
t26 = pkin(2) * t59 + t57 * t76;
t22 = t54 * t123 + t52 * t75;
t18 = pkin(3) * t124 + t29 * t53;
t17 = pkin(3) * t119 + t151 * t53;
t16 = pkin(3) * t120 + t152 * t53;
t15 = pkin(3) * t122 + t153 * t53;
t14 = t74 * (t34 * t69 - t116) * t53 + t34 * t124;
t9 = -t32 * t52 + t86 * t54;
t8 = -t31 * t52 + t87 * t54;
t7 = -t30 * t52 + t88 * t54;
t6 = -t26 * t52 + t89 * t54;
t5 = 0.1e1 / (pkin(2) * t124 + t129 * t139 + t18 * t74);
t4 = 0.1e1 / (pkin(2) * t119 + t53 * t109 + t17 * t72);
t3 = 0.1e1 / (pkin(2) * t120 + t53 * t110 + t16 * t70);
t2 = 0.1e1 / (pkin(2) * t122 + t53 * t111 + t15 * t58);
t1 = -t54 * t34 * t130 + ((t34 * t128 + t33) * t69 + (t34 * t52 - t102) * t75) * t74;
t10 = [((t39 * t129 + t22 * t43) * t139 + (t18 * t39 - t43 * t9) * t74 + (-t43 * t134 + t39 * t55) * t143) * t5, (-(-t43 * t129 + t22 * t39) * t139 + (t18 * t43 + t39 * t9) * t74 + (t39 * t134 + t43 * t55) * t143) * t5, (-(t52 * t123 - t54 * t75) * t139 + (t32 * t54 + t86 * t52) * t74 + t68 * t115) * t5, (-(t1 * t43 + t14 * t39) * (t46 * t85 + t47 * t81) - (-t46 * t81 + t47 * t85) * (t1 * t39 - t14 * t43)) / t14; ((t38 * t131 + t21 * t42) * t140 + (t17 * t38 - t42 * t8) * t72 + (-t42 * t134 + t38 * t55) * t144) * t4, (-(-t42 * t131 + t21 * t38) * t140 + (t17 * t42 + t38 * t8) * t72 + (t38 * t134 + t42 * t55) * t144) * t4, (-(t52 * t125 - t54 * t73) * t140 + (t31 * t54 + t87 * t52) * t72 + t66 * t115) * t4, (-((t150 * t38 - t42 * t90) * t53 + t42 * t112 + (pkin(3) * t106 + t93 * t42) * t72 + pkin(2) * t106) * (t46 * t84 + t47 * t80) - (-t46 * t80 + t47 * t84) * ((-t150 * t42 - t38 * t90) * t53 + t38 * t112 + (-pkin(3) * t103 + t93 * t38) * t72 - pkin(2) * t103)) / (t72 * (t67 * t135 + t151) * t53 + t99 * t119); ((t37 * t132 + t20 * t41) * t141 + (t16 * t37 - t41 * t7) * t70 + (-t41 * t134 + t37 * t55) * t145) * t3, (-(-t41 * t132 + t20 * t37) * t141 + (t16 * t41 + t37 * t7) * t70 + (t37 * t134 + t41 * t55) * t145) * t3, (-(t52 * t126 - t54 * t71) * t141 + (t30 * t54 + t88 * t52) * t70 + t64 * t115) * t3, (-((t149 * t37 - t41 * t91) * t53 + t41 * t113 + (pkin(3) * t107 + t94 * t41) * t70 + pkin(2) * t107) * (t46 * t83 + t47 * t79) - (-t46 * t79 + t47 * t83) * ((-t149 * t41 - t37 * t91) * t53 + t37 * t113 + (-pkin(3) * t104 + t94 * t37) * t70 - pkin(2) * t104)) / (t70 * (t65 * t136 + t152) * t53 + t100 * t120); ((t36 * t133 + t19 * t40) * t142 + (t15 * t36 - t40 * t6) * t58 + (-t40 * t134 + t36 * t55) * t146) * t2, (-(-t40 * t133 + t19 * t36) * t142 + (t15 * t40 + t36 * t6) * t58 + (t36 * t134 + t40 * t55) * t146) * t2, (-(t52 * t127 - t54 * t59) * t142 + (t26 * t54 + t89 * t52) * t58 + t56 * t115) * t2, (-((t148 * t36 - t40 * t92) * t53 + t40 * t114 + (pkin(3) * t108 + t95 * t40) * t58 + pkin(2) * t108) * (t46 * t82 + t47 * t78) - (-t46 * t78 + t47 * t82) * ((-t148 * t40 - t36 * t92) * t53 + t36 * t114 + (-pkin(3) * t105 + t95 * t36) * t58 - pkin(2) * t105)) / (t58 * (t57 * t137 + t153) * t53 + t101 * t122);];
Jinv  = t10;
