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
% Datum: 2020-08-07 03:32
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRRRR10V2G3P3A2_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(8,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR10V2G3P3A2_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR10V2G3P3A2_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3RRRRR10V2G3P3A2_Jinv: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR10V2G3P3A2_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR10V2G3P3A2_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:30:56
% EndTime: 2020-08-07 03:30:58
% DurationCPUTime: 1.62s
% Computational Cost: add. (480->193), mult. (1032->363), div. (9->3), fcn. (903->26), ass. (0->143)
t48 = sin(qJ(3,3));
t105 = t48 * pkin(6) + pkin(3);
t57 = cos(qJ(3,3));
t37 = t57 ^ 2;
t149 = pkin(3) * t37;
t58 = cos(qJ(2,3));
t66 = pkin(8) + pkin(7);
t19 = t66 * t58;
t49 = sin(qJ(2,3));
t50 = sin(qJ(1,3));
t174 = ((pkin(2) * t57 - t105 + 0.2e1 * t149) * t49 - t57 * t19) * t50;
t51 = sin(qJ(3,2));
t104 = t51 * pkin(6) + pkin(3);
t60 = cos(qJ(3,2));
t39 = t60 ^ 2;
t148 = pkin(3) * t39;
t61 = cos(qJ(2,2));
t20 = t66 * t61;
t52 = sin(qJ(2,2));
t53 = sin(qJ(1,2));
t173 = ((pkin(2) * t60 - t104 + 0.2e1 * t148) * t52 - t60 * t20) * t53;
t54 = sin(qJ(3,1));
t103 = t54 * pkin(6) + pkin(3);
t63 = cos(qJ(3,1));
t41 = t63 ^ 2;
t142 = t41 * pkin(3);
t64 = cos(qJ(2,1));
t21 = t66 * t64;
t55 = sin(qJ(2,1));
t56 = sin(qJ(1,1));
t172 = ((pkin(2) * t63 - t103 + 0.2e1 * t142) * t55 - t63 * t21) * t56;
t33 = pkin(2) * t58;
t13 = t33 + pkin(1);
t38 = t58 ^ 2;
t171 = t13 * t49 + (-t38 + 0.1e1) * t66;
t34 = pkin(2) * t61;
t14 = t34 + pkin(1);
t40 = t61 ^ 2;
t170 = t14 * t52 + (-t40 + 0.1e1) * t66;
t35 = pkin(2) * t64;
t15 = t35 + pkin(1);
t42 = t64 ^ 2;
t169 = t15 * t55 + (-t42 + 0.1e1) * t66;
t141 = t42 * pkin(2);
t76 = t55 * t21 - pkin(2) + t141;
t168 = t76 * t54;
t143 = t40 * pkin(2);
t77 = t52 * t20 - pkin(2) + t143;
t167 = t77 * t51;
t144 = t38 * pkin(2);
t78 = t49 * t19 - pkin(2) + t144;
t166 = t78 * t48;
t88 = t104 * t52;
t89 = t105 * t49;
t90 = t55 * t103;
t165 = -pkin(2) * t49 + t19;
t164 = -pkin(2) * t52 + t20;
t163 = -pkin(2) * t55 + t21;
t43 = sin(pkin(4));
t65 = cos(qJ(1,1));
t18 = t66 * t55;
t113 = t18 + t35;
t9 = pkin(1) + t113;
t145 = pkin(3) * t54;
t91 = (t42 - 0.2e1) * t145 - pkin(6);
t156 = t91 * t56 * t43 - t9 * t65;
t62 = cos(qJ(1,2));
t17 = t66 * t52;
t114 = t17 + t34;
t8 = pkin(1) + t114;
t146 = pkin(3) * t51;
t92 = (t40 - 0.2e1) * t146 - pkin(6);
t155 = t92 * t53 * t43 - t8 * t62;
t59 = cos(qJ(1,3));
t16 = t66 * t49;
t115 = t16 + t33;
t7 = pkin(1) + t115;
t147 = pkin(3) * t48;
t93 = (t38 - 0.2e1) * t147 - pkin(6);
t154 = t93 * t50 * t43 - t7 * t59;
t44 = cos(pkin(4));
t153 = t44 / 0.2e1;
t137 = (t44 + 0.1e1) * (t44 - 0.1e1);
t136 = t43 * t49;
t134 = t43 * t52;
t132 = t43 * t55;
t130 = t44 * t57;
t129 = t44 * t59;
t128 = t44 * t60;
t127 = t44 * t62;
t126 = t44 * t63;
t125 = t44 * t65;
t124 = t49 * t50;
t123 = t52 * t53;
t122 = t55 * t56;
t121 = t58 * t59;
t120 = t59 * t49;
t119 = t61 * t62;
t118 = t62 * t52;
t117 = t64 * t65;
t116 = t65 * t55;
t36 = t44 ^ 2;
t112 = t36 - 0.1e1 / 0.2e1;
t111 = pkin(3) * t124;
t110 = pkin(3) * t123;
t109 = pkin(3) * t122;
t45 = legFrame(3,2);
t27 = sin(t45);
t102 = t27 * t111;
t30 = cos(t45);
t101 = t30 * t111;
t46 = legFrame(2,2);
t28 = sin(t46);
t100 = t28 * t110;
t31 = cos(t46);
t99 = t31 * t110;
t47 = legFrame(1,2);
t29 = sin(t47);
t98 = t29 * t109;
t32 = cos(t47);
t97 = t32 * t109;
t96 = t43 * t58 * t120;
t95 = t43 * t61 * t118;
t94 = t43 * t64 * t116;
t87 = (t16 + pkin(1)) * t58 + t144;
t86 = (t17 + pkin(1)) * t61 + t143;
t85 = (t18 + pkin(1)) * t64 + t141;
t84 = t165 * t43;
t83 = t164 * t43;
t82 = t163 * t43;
t81 = t171 * t43 * t59 + pkin(6) * t124;
t80 = t43 * t170 * t62 + pkin(6) * t123;
t79 = t43 * t169 * t65 + pkin(6) * t122;
t72 = t50 * t166;
t71 = t53 * t167;
t70 = t56 * t168;
t69 = -t93 * t57 - t166;
t68 = -t92 * t60 - t167;
t67 = -t91 * t63 - t168;
t3 = 0.1e1 / ((-t64 * pkin(6) * t142 + (-t113 * pkin(6) + t15 * t145) * t63 + t54 * pkin(2) * t9) * t43 + (t21 + (-pkin(3) * t63 - pkin(2)) * t55) * pkin(1) * t126);
t2 = 0.1e1 / ((-t61 * pkin(6) * t148 + (-t114 * pkin(6) + t14 * t146) * t60 + t51 * pkin(2) * t8) * t43 + (t20 + (-pkin(3) * t60 - pkin(2)) * t52) * pkin(1) * t128);
t1 = 0.1e1 / ((-t58 * pkin(6) * t149 + (-t115 * pkin(6) + t13 * t147) * t57 + t48 * pkin(2) * t7) * t43 + (t19 + (-pkin(3) * t57 - pkin(2)) * t49) * pkin(1) * t130);
t4 = [((-t32 * t172 + t67 * t29) * t36 + ((t32 * t117 + 0.2e1 * t29 * t132) * t142 + (-t156 * t32 - t29 * t82) * t63 - (t29 * t90 + t32 * t70) * t43) * t44 + t41 * t97 + (-(-t29 * t42 + t32 * t94 + t29) * t145 - t29 * pkin(6)) * t63 + (t85 * t29 - t79 * t32) * t54 - t97) * t3, ((t29 * t172 + t67 * t32) * t36 + (-(t29 * t117 - 0.2e1 * t32 * t132) * t142 + (t156 * t29 - t32 * t82) * t63 + (t29 * t70 - t32 * t90) * t43) * t44 - t41 * t98 + ((t29 * t94 + t32 * t42 - t32) * t145 - t32 * pkin(6)) * t63 + (t79 * t29 + t85 * t32) * t54 + t98) * t3, (((t64 * t54 * t109 - t91 * t125) * t63 + (-t76 * t125 + t56 * t169) * t54) * t43 - 0.2e1 * (t56 * t64 * t153 + t112 * t116) * t142 - (-t163 * t125 + t56 * t9) * t126 + t65 * t90 * t137) * t3; ((-t31 * t173 + t68 * t28) * t36 + ((t31 * t119 + 0.2e1 * t28 * t134) * t148 + (-t155 * t31 - t28 * t83) * t60 - (t28 * t88 + t31 * t71) * t43) * t44 + t39 * t99 + (-(-t28 * t40 + t31 * t95 + t28) * t146 - t28 * pkin(6)) * t60 + (t86 * t28 - t80 * t31) * t51 - t99) * t2, ((t28 * t173 + t68 * t31) * t36 + (-(t28 * t119 - 0.2e1 * t31 * t134) * t148 + (t155 * t28 - t31 * t83) * t60 + (t28 * t71 - t31 * t88) * t43) * t44 - t39 * t100 + ((t28 * t95 + t31 * t40 - t31) * t146 - t31 * pkin(6)) * t60 + (t80 * t28 + t86 * t31) * t51 + t100) * t2, (((t61 * t51 * t110 - t92 * t127) * t60 + (-t77 * t127 + t53 * t170) * t51) * t43 - 0.2e1 * (t53 * t61 * t153 + t112 * t118) * t148 - (-t164 * t127 + t53 * t8) * t128 + t62 * t88 * t137) * t2; ((-t30 * t174 + t69 * t27) * t36 + ((t30 * t121 + 0.2e1 * t27 * t136) * t149 + (-t154 * t30 - t27 * t84) * t57 - (t27 * t89 + t30 * t72) * t43) * t44 + t37 * t101 + (-(-t27 * t38 + t30 * t96 + t27) * t147 - t27 * pkin(6)) * t57 + (t87 * t27 - t81 * t30) * t48 - t101) * t1, ((t27 * t174 + t69 * t30) * t36 + (-(t27 * t121 - 0.2e1 * t30 * t136) * t149 + (t154 * t27 - t30 * t84) * t57 + (t27 * t72 - t30 * t89) * t43) * t44 - t37 * t102 + ((t27 * t96 + t30 * t38 - t30) * t147 - t30 * pkin(6)) * t57 + (t81 * t27 + t87 * t30) * t48 + t102) * t1, (((t58 * t48 * t111 - t93 * t129) * t57 + (-t78 * t129 + t50 * t171) * t48) * t43 - 0.2e1 * (t50 * t58 * t153 + t112 * t120) * t149 - (-t165 * t129 + t50 * t7) * t130 + t59 * t89 * t137) * t1;];
Jinv  = t4;
