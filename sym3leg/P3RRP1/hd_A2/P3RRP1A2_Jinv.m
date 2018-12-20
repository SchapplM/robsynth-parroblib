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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
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
% Datum: 2018-12-20 18:11
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Jinv = P3RRP1A2_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(4,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRP1A2_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRP1A2_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRP1A2_Jinv: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRP1A2_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRP1A2_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 18:11:29
% EndTime: 2018-12-20 18:11:31
% DurationCPUTime: 1.72s
% Computational Cost: add. (702->244), mult. (1710->416), div. (9->6), fcn. (603->20), ass. (0->194)
t134 = xP(3);
t100 = sin(t134);
t101 = cos(t134);
t140 = koppelP(1,2);
t143 = koppelP(1,1);
t203 = (pkin(2) * qJ(3,1));
t172 = -2 * t203;
t137 = (qJ(3,1) ^ 2);
t144 = (pkin(2) ^ 2);
t96 = (-t137 + t144);
t47 = t140 * t172 + t96 * t143;
t104 = qJ(3,1) * t143;
t210 = 2 * pkin(2);
t48 = t104 * t210 + t96 * t140;
t158 = t100 * t48 - t101 * t47;
t24 = t100 * t47 + t101 * t48;
t127 = legFrame(1,3);
t90 = sin(t127);
t93 = cos(t127);
t216 = t158 * t93 - t90 * t24;
t139 = koppelP(2,2);
t142 = koppelP(2,1);
t202 = pkin(2) * qJ(3,2);
t171 = -2 * t202;
t136 = qJ(3,2) ^ 2;
t95 = -t136 + t144;
t45 = t139 * t171 + t95 * t142;
t103 = qJ(3,2) * t142;
t46 = t103 * t210 + t95 * t139;
t159 = t100 * t46 - t101 * t45;
t22 = t100 * t45 + t101 * t46;
t126 = legFrame(2,3);
t89 = sin(t126);
t92 = cos(t126);
t215 = t159 * t92 - t89 * t22;
t138 = koppelP(3,2);
t141 = koppelP(3,1);
t201 = pkin(2) * qJ(3,3);
t170 = -2 * t201;
t135 = qJ(3,3) ^ 2;
t94 = -t135 + t144;
t43 = t138 * t170 + t94 * t141;
t102 = qJ(3,3) * t141;
t44 = t102 * t210 + t94 * t138;
t160 = t100 * t44 - t101 * t43;
t20 = t100 * t43 + t101 * t44;
t125 = legFrame(3,3);
t88 = sin(t125);
t91 = cos(t125);
t214 = t160 * t91 - t88 * t20;
t131 = cos(qJ(1,3));
t114 = t131 ^ 2;
t145 = (pkin(1) ^ 2);
t120 = 2 + t145;
t173 = 2 * pkin(1);
t128 = sin(qJ(1,3));
t181 = t128 * t131;
t183 = qJ(3,3) * t128;
t117 = qJ(1,3) + qJ(2,3);
t82 = sin(t117);
t209 = pkin(1) * t82;
t79 = 1 + t94;
t85 = cos(t117);
t213 = -(t131 * t82 * qJ(3,3) - ((t79 * t181 + (0.2e1 * t114 - 0.1e1) * t201) * t209 + t183) * t85) * t173 - (t135 * t120);
t132 = cos(qJ(1,2));
t115 = t132 ^ 2;
t129 = sin(qJ(1,2));
t180 = t129 * t132;
t186 = qJ(3,2) * t129;
t118 = qJ(1,2) + qJ(2,2);
t83 = sin(t118);
t208 = pkin(1) * t83;
t80 = 1 + t95;
t86 = cos(t118);
t212 = -(t132 * t83 * qJ(3,2) - ((t80 * t180 + (0.2e1 * t115 - 0.1e1) * t202) * t208 + t186) * t86) * t173 - (t136 * t120);
t133 = cos(qJ(1,1));
t116 = t133 ^ 2;
t130 = sin(qJ(1,1));
t179 = t130 * t133;
t189 = qJ(3,1) * t130;
t119 = qJ(1,1) + qJ(2,1);
t84 = sin(t119);
t207 = pkin(1) * t84;
t81 = 1 + t96;
t87 = cos(t119);
t211 = -(t133 * t84 * qJ(3,1) - ((t81 * t179 + (0.2e1 * t116 - 0.1e1) * t203) * t207 + t189) * t87) * t173 - (t137 * t120);
t206 = pkin(2) * t91;
t205 = pkin(2) * t92;
t204 = pkin(2) * t93;
t73 = t88 * pkin(2);
t74 = t89 * pkin(2);
t75 = t90 * pkin(2);
t200 = pkin(2) * t145;
t199 = t82 * t85;
t198 = t83 * t86;
t197 = t84 * t87;
t193 = qJ(3,1) * t93;
t192 = qJ(3,2) * t92;
t191 = qJ(3,3) * t91;
t70 = t88 * qJ(3,3);
t71 = t89 * qJ(3,2);
t72 = t90 * qJ(3,1);
t121 = 1 + t145;
t190 = qJ(3,1) * t121;
t188 = qJ(3,1) * t140;
t187 = qJ(3,2) * t121;
t185 = qJ(3,2) * t139;
t184 = qJ(3,3) * t121;
t182 = qJ(3,3) * t138;
t175 = -0.1e1 / 0.2e1 - t144 / 0.2e1;
t174 = 0.1e1 / 0.4e1 + t144 / 0.4e1;
t169 = pkin(2) * t70;
t168 = pkin(2) * t71;
t167 = pkin(2) * t72;
t166 = t135 / 0.2e1 + t175;
t165 = t136 / 0.2e1 + t175;
t164 = t137 / 0.2e1 + t175;
t163 = t179 * t203;
t162 = t180 * t202;
t161 = t181 * t201;
t37 = t70 - t206;
t38 = t73 + t191;
t154 = -t114 * t37 - t38 * t181;
t153 = t114 * t38 - t37 * t181;
t39 = t71 - t205;
t40 = t74 + t192;
t152 = -t115 * t39 - t40 * t180;
t151 = t115 * t40 - t39 * t180;
t41 = t72 - t204;
t42 = t75 + t193;
t150 = -t116 * t41 - t42 * t179;
t149 = t116 * t42 - t41 * t179;
t148 = -t81 * t116 + 0.2e1 * t163;
t147 = -t80 * t115 + 0.2e1 * t162;
t146 = -t79 * t114 + 0.2e1 * t161;
t123 = t144 / 0.2e1;
t110 = pkin(2) * t143;
t109 = pkin(2) * t142;
t108 = pkin(2) * t141;
t107 = pkin(2) * t140;
t106 = pkin(2) * t139;
t105 = pkin(2) * t138;
t78 = t87 ^ 2;
t77 = t86 ^ 2;
t76 = t85 ^ 2;
t69 = -t104 + t107;
t68 = t104 + t107;
t67 = t110 - t188;
t66 = t110 + t188;
t65 = -t103 + t106;
t64 = t103 + t106;
t63 = t109 - t185;
t62 = t109 + t185;
t61 = -t102 + t105;
t60 = t102 + t105;
t59 = t108 - t182;
t58 = t108 + t182;
t57 = t140 * t200 + t104;
t56 = t143 * t200 - t188;
t55 = t139 * t200 + t103;
t54 = t142 * t200 - t185;
t53 = t138 * t200 + t102;
t52 = t141 * t200 - t182;
t51 = -pkin(2) * t133 + t189;
t50 = -pkin(2) * t132 + t186;
t49 = -pkin(2) * t131 + t183;
t36 = t93 * t172 + t90 * t96;
t35 = t92 * t171 + t89 * t95;
t34 = t91 * t170 + t88 * t94;
t33 = t100 * t66 + t101 * t69;
t32 = -t100 * t69 + t101 * t66;
t31 = t100 * t62 + t101 * t65;
t30 = -t100 * t65 + t101 * t62;
t29 = t100 * t58 + t101 * t61;
t28 = -t100 * t61 + t101 * t58;
t18 = (t93 * t96 + 0.2e1 * t167) * t133 + t36 * t130;
t17 = (t92 * t95 + 0.2e1 * t168) * t132 + t35 * t129;
t16 = (t91 * t94 + 0.2e1 * t169) * t131 + t34 * t128;
t15 = t36 * t133 - 0.2e1 * t130 * ((t123 - t137 / 0.2e1) * t93 + t167);
t14 = t35 * t132 - 0.2e1 * t129 * ((t123 - t136 / 0.2e1) * t92 + t168);
t13 = t34 * t131 - 0.2e1 * t128 * ((t123 - t135 / 0.2e1) * t91 + t169);
t12 = -t32 * t90 + t33 * t93;
t11 = t32 * t93 + t33 * t90;
t10 = -t30 * t89 + t31 * t92;
t9 = t30 * t92 + t31 * t89;
t8 = -t28 * t88 + t29 * t91;
t7 = t28 * t91 + t29 * t88;
t6 = t158 * t90 + t24 * t93;
t5 = t159 * t89 + t22 * t92;
t4 = t160 * t88 + t20 * t91;
t3 = 0.1e1 / ((0.2e1 * (-t148 + t164) * t78 + t148) * t145 + t211);
t2 = 0.1e1 / ((0.2e1 * (-t147 + t165) * t77 + t147) * t145 + t212);
t1 = 0.1e1 / ((0.2e1 * (-t146 + t166) * t76 + t146) * t145 + t213);
t19 = [((t121 * t90 * t84 - t93 * t87) * qJ(3,1) + (t15 * t78 + t18 * t197 + t51 * (t75 - t193)) * pkin(1) + ((-t149 + t75) * t87 + t150 * t84) * t145) * t3 (-t93 * t84 * t190 - t72 * t87 + (-t18 * t78 + t15 * t197 - t51 * (t72 + t204)) * pkin(1) + ((t150 - t204) * t87 + t149 * t84) * t145) * t3 (((-t130 * t216 + t6 * t133) * t207 + (-t100 * t57 + t101 * t56) * t93 + t90 * (t100 * t56 + t101 * t57) + (-t11 * t116 - t12 * t179) * t145) * t87 + (((-t100 * t140 + t101 * t143) * t93 + t90 * (t100 * t143 + t101 * t140)) * t190 + (-t11 * t179 + t12 * t116) * t145) * t84 + (-(t130 * t6 + t216 * t133) * t78 + ((-t100 * t68 + t101 * t67) * t93 + (t100 * t67 + t101 * t68) * t90) * t51) * pkin(1)) / ((0.4e1 * (t164 * t116 + t163 - t137 / 0.4e1 + t174) * t78 - t148) * t145 - t211); ((t121 * t89 * t83 - t92 * t86) * qJ(3,2) + (t14 * t77 + t17 * t198 + t50 * (t74 - t192)) * pkin(1) + ((-t151 + t74) * t86 + t152 * t83) * t145) * t2 (-t92 * t83 * t187 - t71 * t86 + (-t17 * t77 + t14 * t198 - t50 * (t71 + t205)) * pkin(1) + ((t152 - t205) * t86 + t151 * t83) * t145) * t2 (((-t129 * t215 + t5 * t132) * t208 + (-t100 * t55 + t101 * t54) * t92 + t89 * (t100 * t54 + t101 * t55) + (-t10 * t180 - t9 * t115) * t145) * t86 + (((-t100 * t139 + t101 * t142) * t92 + t89 * (t100 * t142 + t101 * t139)) * t187 + (t10 * t115 - t9 * t180) * t145) * t83 + (-(t129 * t5 + t215 * t132) * t77 + ((-t100 * t64 + t101 * t63) * t92 + (t100 * t63 + t101 * t64) * t89) * t50) * pkin(1)) / ((0.4e1 * (t165 * t115 + t162 - t136 / 0.4e1 + t174) * t77 - t147) * t145 - t212); ((t121 * t88 * t82 - t91 * t85) * qJ(3,3) + (t13 * t76 + t16 * t199 + t49 * (t73 - t191)) * pkin(1) + ((-t153 + t73) * t85 + t154 * t82) * t145) * t1 (-t91 * t82 * t184 - t70 * t85 + (-t16 * t76 + t13 * t199 - t49 * (t70 + t206)) * pkin(1) + ((t154 - t206) * t85 + t153 * t82) * t145) * t1 (((-t128 * t214 + t4 * t131) * t209 + (-t100 * t53 + t101 * t52) * t91 + t88 * (t100 * t52 + t101 * t53) + (-t7 * t114 - t8 * t181) * t145) * t85 + (((-t100 * t138 + t101 * t141) * t91 + t88 * (t100 * t141 + t101 * t138)) * t184 + (t8 * t114 - t7 * t181) * t145) * t82 + (-(t128 * t4 + t214 * t131) * t76 + ((-t100 * t60 + t101 * t59) * t91 + (t100 * t59 + t101 * t60) * t88) * t49) * pkin(1)) / ((0.4e1 * (t166 * t114 + t161 - t135 / 0.4e1 + t174) * t76 - t146) * t145 - t213);];
Jinv  = t19;
