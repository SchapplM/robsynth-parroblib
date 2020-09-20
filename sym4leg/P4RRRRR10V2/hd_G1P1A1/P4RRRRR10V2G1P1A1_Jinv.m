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
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
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
% Datum: 2020-08-07 17:23
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P4RRRRR10V2G1P1A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(8,1),zeros(4,3),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4RRRRR10V2G1P1A1_Jinv: qJ has to be [3x4] (double)');
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4RRRRR10V2G1P1A1_Jinv: xP has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4RRRRR10V2G1P1A1_Jinv: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4RRRRR10V2G1P1A1_Jinv: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4RRRRR10V2G1P1A1_Jinv: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 17:20:37
% EndTime: 2020-08-07 17:20:38
% DurationCPUTime: 1.25s
% Computational Cost: add. (680->242), mult. (1396->485), div. (16->8), fcn. (1168->36), ass. (0->187)
t86 = cos(pkin(4));
t194 = pkin(1) * t86;
t85 = sin(pkin(4));
t193 = pkin(2) * t85;
t192 = pkin(2) * t86;
t95 = cos(qJ(2,4));
t191 = pkin(2) * t95;
t94 = cos(qJ(3,4));
t81 = t94 ^ 2;
t190 = pkin(3) * t81;
t106 = cos(qJ(3,3));
t82 = t106 ^ 2;
t189 = pkin(3) * t82;
t109 = cos(qJ(3,2));
t83 = t109 ^ 2;
t188 = pkin(3) * t83;
t112 = cos(qJ(3,1));
t84 = t112 ^ 2;
t187 = pkin(3) * t84;
t186 = pkin(3) * t85;
t91 = sin(qJ(3,4));
t185 = pkin(3) * t91;
t97 = sin(qJ(3,3));
t184 = pkin(3) * t97;
t183 = pkin(6) * t85;
t107 = cos(qJ(2,3));
t182 = pkin(2) * t107;
t110 = cos(qJ(2,2));
t181 = pkin(2) * t110;
t113 = cos(qJ(2,1));
t180 = pkin(2) * t113;
t100 = sin(qJ(3,2));
t179 = pkin(3) * t100;
t103 = sin(qJ(3,1));
t178 = pkin(3) * t103;
t177 = t85 * t91;
t176 = t85 * t97;
t175 = t86 * t91;
t92 = sin(qJ(2,4));
t174 = t86 * t92;
t173 = t86 * t97;
t98 = sin(qJ(2,3));
t172 = t86 * t98;
t171 = t100 * t85;
t170 = t100 * t86;
t101 = sin(qJ(2,2));
t169 = t101 * t86;
t168 = t103 * t85;
t167 = t103 * t86;
t104 = sin(qJ(2,1));
t166 = t104 * t86;
t115 = pkin(8) + pkin(7);
t165 = t115 * t86;
t164 = t115 * t92;
t163 = t115 * t98;
t117 = koppelP(4,2);
t162 = t117 * t86;
t118 = koppelP(3,2);
t161 = t118 * t86;
t119 = koppelP(2,2);
t160 = t119 * t86;
t120 = koppelP(1,2);
t159 = t120 * t86;
t121 = koppelP(4,1);
t158 = t121 * t86;
t122 = koppelP(3,1);
t157 = t122 * t86;
t123 = koppelP(2,1);
t156 = t123 * t86;
t124 = koppelP(1,1);
t155 = t124 * t86;
t154 = t115 * t101;
t153 = t115 * t104;
t152 = t115 * t117;
t151 = t115 * t118;
t150 = t115 * t119;
t149 = t115 * t120;
t148 = t115 * t121;
t147 = t115 * t122;
t146 = t115 * t123;
t145 = t115 * t124;
t144 = pkin(3) * t177;
t143 = pkin(3) * t176;
t142 = pkin(2) * t177;
t141 = pkin(2) * t176;
t140 = pkin(3) * t171;
t139 = pkin(3) * t168;
t138 = pkin(2) * t171;
t137 = pkin(2) * t168;
t87 = legFrame(4,3);
t71 = sin(t87);
t75 = cos(t87);
t93 = sin(qJ(1,4));
t96 = cos(qJ(1,4));
t33 = t71 * t96 + t75 * t93;
t136 = t33 * t177;
t34 = -t71 * t93 + t75 * t96;
t135 = t34 * t177;
t108 = cos(qJ(1,3));
t88 = legFrame(3,3);
t72 = sin(t88);
t76 = cos(t88);
t99 = sin(qJ(1,3));
t35 = t108 * t72 + t76 * t99;
t134 = t35 * t176;
t36 = t108 * t76 - t72 * t99;
t133 = t36 * t176;
t102 = sin(qJ(1,2));
t111 = cos(qJ(1,2));
t89 = legFrame(2,3);
t73 = sin(t89);
t77 = cos(t89);
t37 = t102 * t77 + t111 * t73;
t132 = t37 * t171;
t38 = -t102 * t73 + t111 * t77;
t131 = t38 * t171;
t105 = sin(qJ(1,1));
t114 = cos(qJ(1,1));
t90 = legFrame(1,3);
t74 = sin(t90);
t78 = cos(t90);
t39 = t105 * t78 + t114 * t74;
t130 = t39 * t168;
t40 = -t105 * t74 + t114 * t78;
t129 = t40 * t168;
t128 = -(pkin(1) * t174 + t95 * t183) * t190 + (pkin(1) + t164 + t191) * t142;
t127 = -(pkin(1) * t172 + t107 * t183) * t189 + (pkin(1) + t163 + t182) * t141;
t126 = -(pkin(1) * t169 + t110 * t183) * t188 + (pkin(1) + t154 + t181) * t138;
t125 = -(pkin(1) * t166 + t113 * t183) * t187 + (pkin(1) + t153 + t180) * t137;
t116 = xP(4);
t80 = cos(t116);
t79 = sin(t116);
t70 = -pkin(6) + t178;
t69 = -pkin(6) + t179;
t68 = -pkin(6) + t184;
t67 = -pkin(6) + t185;
t66 = pkin(1) * t165;
t65 = -pkin(2) * t104 + t113 * t115;
t64 = -pkin(2) * t101 + t110 * t115;
t63 = -pkin(2) * t98 + t107 * t115;
t62 = -pkin(2) * t92 + t115 * t95;
t61 = t105 * t120 + t114 * t124;
t60 = t105 * t124 - t114 * t120;
t59 = t102 * t119 + t111 * t123;
t58 = t102 * t123 - t111 * t119;
t57 = t108 * t122 + t118 * t99;
t56 = -t108 * t118 + t122 * t99;
t55 = t117 * t93 + t121 * t96;
t54 = -t117 * t96 + t121 * t93;
t53 = -pkin(1) * t192 - t115 * t183;
t52 = t104 * t155 - t120 * t113;
t51 = t101 * t156 - t119 * t110;
t50 = -t118 * t107 + t98 * t157;
t49 = t104 * t159 + t124 * t113;
t48 = t101 * t160 + t123 * t110;
t47 = t122 * t107 + t98 * t161;
t46 = -t117 * t95 + t92 * t158;
t45 = t121 * t95 + t92 * t162;
t28 = t105 * t49 + t114 * t52;
t27 = t102 * t48 + t111 * t51;
t26 = t108 * t50 + t47 * t99;
t25 = -t105 * t52 + t114 * t49;
t24 = -t102 * t51 + t111 * t48;
t23 = t108 * t47 - t50 * t99;
t22 = t45 * t93 + t46 * t96;
t21 = t45 * t96 - t46 * t93;
t20 = (pkin(2) * t120 + t86 * t145) * t113 + (-pkin(2) * t155 + t149) * t104 + t124 * t139;
t19 = (pkin(2) * t119 + t86 * t146) * t110 + (-pkin(2) * t156 + t150) * t101 + t123 * t140;
t18 = (pkin(2) * t118 + t86 * t147) * t107 + (-pkin(2) * t157 + t151) * t98 + t122 * t143;
t17 = (pkin(2) * t124 - t86 * t149) * t113 + (pkin(2) * t159 + t145) * t104 - t120 * t139;
t16 = (pkin(2) * t123 - t86 * t150) * t110 + (pkin(2) * t160 + t146) * t101 - t119 * t140;
t15 = (pkin(2) * t122 - t86 * t151) * t107 + (pkin(2) * t161 + t147) * t98 - t118 * t143;
t14 = (pkin(2) * t117 + t86 * t148) * t95 + (-pkin(2) * t158 + t152) * t92 + t121 * t144;
t13 = (pkin(2) * t121 - t86 * t152) * t95 + (pkin(2) * t162 + t148) * t92 - t117 * t144;
t12 = 0.1e1 / (((pkin(1) * t178 - pkin(6) * t153 + t70 * t180) * t85 + t65 * t194) * t112 + t125);
t11 = 0.1e1 / (((pkin(1) * t179 - pkin(6) * t154 + t69 * t181) * t85 + t64 * t194) * t109 + t126);
t10 = 0.1e1 / (((pkin(1) * t184 - pkin(6) * t163 + t68 * t182) * t85 + t63 * t194) * t106 + t127);
t9 = 0.1e1 / (((pkin(1) * t185 - pkin(6) * t164 + t67 * t191) * t85 + t62 * t194) * t94 + t128);
t8 = -t105 * t17 + t114 * t20;
t7 = -t102 * t16 + t111 * t19;
t6 = t108 * t18 - t15 * t99;
t5 = t105 * t20 + t114 * t17;
t4 = t102 * t19 + t111 * t16;
t3 = t108 * t15 + t18 * t99;
t2 = -t13 * t93 + t14 * t96;
t1 = t13 * t96 + t14 * t93;
t29 = [(-(t113 * t40 - t39 * t166) * t187 + (-pkin(3) * t130 + (-pkin(2) * t40 - t39 * t165) * t113 + t104 * (-t115 * t40 + t39 * t192)) * t112 - pkin(2) * t130) * t12, (-(t113 * t39 + t40 * t166) * t187 + (pkin(3) * t129 + (-pkin(2) * t39 + t40 * t165) * t113 - (t115 * t39 + t40 * t192) * t104) * t112 + pkin(2) * t129) * t12, (-t104 * t84 * t186 + (-pkin(3) * t167 + t65 * t85) * t112 - pkin(2) * t167) * t12, (-((-t25 * t79 + t28 * t80) * t78 + t74 * (t25 * t80 + t28 * t79)) * t187 + ((t5 * t79 + t8 * t80) * t78 - t74 * (t5 * t80 - t79 * t8)) * t112 + ((t60 * t79 + t61 * t80) * t78 + t74 * (-t60 * t80 + t61 * t79)) * t137) / (((t70 * t193 + t66) * t113 + t53 * t104 + pkin(1) * t139) * t112 + t125); (-(t110 * t38 - t37 * t169) * t188 + (-pkin(3) * t132 + (-pkin(2) * t38 - t37 * t165) * t110 + t101 * (-t115 * t38 + t37 * t192)) * t109 - pkin(2) * t132) * t11, (-(t110 * t37 + t38 * t169) * t188 + (pkin(3) * t131 + (-pkin(2) * t37 + t38 * t165) * t110 - (t115 * t37 + t38 * t192) * t101) * t109 + pkin(2) * t131) * t11, (-t101 * t83 * t186 + (-pkin(3) * t170 + t64 * t85) * t109 - pkin(2) * t170) * t11, (-((-t24 * t79 + t27 * t80) * t77 + t73 * (t24 * t80 + t27 * t79)) * t188 + ((t4 * t79 + t7 * t80) * t77 - t73 * (t4 * t80 - t7 * t79)) * t109 + ((t58 * t79 + t59 * t80) * t77 + t73 * (-t58 * t80 + t59 * t79)) * t138) / (((t69 * t193 + t66) * t110 + t53 * t101 + pkin(1) * t140) * t109 + t126); (-(t107 * t36 - t35 * t172) * t189 + (-pkin(3) * t134 + (-pkin(2) * t36 - t35 * t165) * t107 + t98 * (-t115 * t36 + t35 * t192)) * t106 - pkin(2) * t134) * t10, (-(t107 * t35 + t36 * t172) * t189 + (pkin(3) * t133 + (-pkin(2) * t35 + t36 * t165) * t107 - (t115 * t35 + t36 * t192) * t98) * t106 + pkin(2) * t133) * t10, (-t98 * t82 * t186 + (-pkin(3) * t173 + t63 * t85) * t106 - pkin(2) * t173) * t10, (-((-t23 * t79 + t26 * t80) * t76 + t72 * (t23 * t80 + t26 * t79)) * t189 + ((t3 * t79 + t6 * t80) * t76 - t72 * (t3 * t80 - t6 * t79)) * t106 + ((t56 * t79 + t57 * t80) * t76 + t72 * (-t56 * t80 + t57 * t79)) * t141) / (((t68 * t193 + t66) * t107 + t53 * t98 + pkin(1) * t143) * t106 + t127); (-(-t33 * t174 + t34 * t95) * t190 + (-pkin(3) * t136 + (-pkin(2) * t34 - t33 * t165) * t95 + t92 * (-t115 * t34 + t33 * t192)) * t94 - pkin(2) * t136) * t9, (-(t34 * t174 + t33 * t95) * t190 + (pkin(3) * t135 + (-pkin(2) * t33 + t34 * t165) * t95 - (t115 * t33 + t34 * t192) * t92) * t94 + pkin(2) * t135) * t9, (-t92 * t81 * t186 + (-pkin(3) * t175 + t62 * t85) * t94 - pkin(2) * t175) * t9, (-((-t21 * t79 + t22 * t80) * t75 + t71 * (t21 * t80 + t22 * t79)) * t190 + ((t1 * t79 + t2 * t80) * t75 - t71 * (t1 * t80 - t2 * t79)) * t94 + ((t54 * t79 + t55 * t80) * t75 + t71 * (-t54 * t80 + t55 * t79)) * t142) / (((t67 * t193 + t66) * t95 + t53 * t92 + pkin(1) * t144) * t94 + t128);];
Jinv  = t29;
