% Analytische Jacobi-Matrix für parallelen Roboter
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorpose und aktiven Gelenkkoordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Input:
% xP [6x1]
%   Generalized platform coordinates
% qJ [3x6]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [6x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,alpha3,alpha4,d4,theta1,theta2,theta3]';
% koppelP [6x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% Jinv [6x6]
%   Analytische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-16 19:51
% Revision: 36f6366a01c4a552c0708fcd8ed3e0fb9da693e2 (2019-05-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P6PPPRRR1V2G1P1A3_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,6),zeros(10,1),zeros(6,3),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [3 6]), ...
  'P6PPPRRR1V2G1P1A3_Jinv: qJ has to be [3x6] (double)');
assert(isreal(xP) && all(size(xP) == [6 1]), ...
  'P6PPPRRR1V2G1P1A3_Jinv: xP has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'P6PPPRRR1V2G1P1A3_Jinv: pkin has to be [10x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [6 3]), ...
  'P6PPPRRR1V2G1P1A3_Jinv: legFrame has to be [6x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6PPPRRR1V2G1P1A3_Jinv: Koppelpunkt has to be [6x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-16 19:50:40
% EndTime: 2019-05-16 19:50:41
% DurationCPUTime: 1.44s
% Computational Cost: add. (528->198), mult. (1242->457), div. (72->2), fcn. (1374->46), ass. (0->215)
t125 = xP(6);
t100 = cos(t125);
t127 = xP(4);
t102 = cos(t127);
t145 = koppelP(1,1);
t165 = t102 * t145;
t139 = koppelP(1,2);
t171 = t102 * t139;
t126 = xP(5);
t101 = cos(t126);
t133 = koppelP(1,3);
t177 = t101 * t133;
t99 = sin(t127);
t189 = t145 * t99;
t195 = t139 * t99;
t97 = sin(t125);
t152 = t100 * t145 - t139 * t97;
t98 = sin(t126);
t30 = t152 * t101 + t98 * t133;
t124 = legFrame(1,2);
t90 = sin(t124);
t223 = t90 * t30;
t118 = legFrame(1,1);
t72 = sin(t118);
t96 = cos(t124);
t240 = (-(t98 * t189 + t171) * t100 + t99 * t177 - t97 * (-t98 * t195 + t165)) * t96 + t72 * t223;
t144 = koppelP(2,1);
t166 = t102 * t144;
t138 = koppelP(2,2);
t172 = t102 * t138;
t132 = koppelP(2,3);
t178 = t101 * t132;
t190 = t144 * t99;
t196 = t138 * t99;
t153 = t100 * t144 - t138 * t97;
t29 = t153 * t101 + t98 * t132;
t123 = legFrame(2,2);
t89 = sin(t123);
t224 = t89 * t29;
t117 = legFrame(2,1);
t71 = sin(t117);
t95 = cos(t123);
t239 = (-(t98 * t190 + t172) * t100 + t99 * t178 - t97 * (-t98 * t196 + t166)) * t95 + t71 * t224;
t143 = koppelP(3,1);
t167 = t102 * t143;
t137 = koppelP(3,2);
t173 = t102 * t137;
t131 = koppelP(3,3);
t179 = t101 * t131;
t191 = t143 * t99;
t197 = t137 * t99;
t154 = t100 * t143 - t137 * t97;
t28 = t154 * t101 + t98 * t131;
t122 = legFrame(3,2);
t88 = sin(t122);
t225 = t88 * t28;
t116 = legFrame(3,1);
t70 = sin(t116);
t94 = cos(t122);
t238 = (-(t98 * t191 + t173) * t100 + t99 * t179 - t97 * (-t98 * t197 + t167)) * t94 + t70 * t225;
t142 = koppelP(4,1);
t168 = t102 * t142;
t136 = koppelP(4,2);
t174 = t102 * t136;
t130 = koppelP(4,3);
t180 = t101 * t130;
t192 = t142 * t99;
t198 = t136 * t99;
t155 = t100 * t142 - t136 * t97;
t27 = t155 * t101 + t98 * t130;
t121 = legFrame(4,2);
t87 = sin(t121);
t226 = t87 * t27;
t115 = legFrame(4,1);
t69 = sin(t115);
t93 = cos(t121);
t237 = (-(t98 * t192 + t174) * t100 + t99 * t180 - t97 * (-t98 * t198 + t168)) * t93 + t69 * t226;
t141 = koppelP(5,1);
t169 = t102 * t141;
t135 = koppelP(5,2);
t175 = t102 * t135;
t129 = koppelP(5,3);
t181 = t101 * t129;
t193 = t141 * t99;
t199 = t135 * t99;
t156 = t100 * t141 - t135 * t97;
t26 = t156 * t101 + t98 * t129;
t120 = legFrame(5,2);
t86 = sin(t120);
t227 = t86 * t26;
t114 = legFrame(5,1);
t68 = sin(t114);
t92 = cos(t120);
t236 = (-(t98 * t193 + t175) * t100 + t99 * t181 - t97 * (-t98 * t199 + t169)) * t92 + t68 * t227;
t140 = koppelP(6,1);
t170 = t102 * t140;
t134 = koppelP(6,2);
t176 = t102 * t134;
t128 = koppelP(6,3);
t182 = t101 * t128;
t194 = t140 * t99;
t200 = t134 * t99;
t157 = t100 * t140 - t134 * t97;
t25 = t157 * t101 + t98 * t128;
t119 = legFrame(6,2);
t85 = sin(t119);
t228 = t85 * t25;
t113 = legFrame(6,1);
t67 = sin(t113);
t91 = cos(t119);
t235 = (-(t98 * t194 + t176) * t100 + t99 * t182 - t97 * (-t98 * t200 + t170)) * t91 + t67 * t228;
t105 = sin(pkin(8));
t106 = cos(pkin(8));
t107 = legFrame(6,3);
t61 = sin(t107);
t73 = cos(t107);
t49 = -t61 * t105 + t106 * t73;
t234 = t49 * t85;
t108 = legFrame(5,3);
t62 = sin(t108);
t74 = cos(t108);
t50 = -t62 * t105 + t106 * t74;
t233 = t50 * t86;
t109 = legFrame(4,3);
t63 = sin(t109);
t75 = cos(t109);
t51 = -t63 * t105 + t106 * t75;
t232 = t51 * t87;
t110 = legFrame(3,3);
t64 = sin(t110);
t76 = cos(t110);
t52 = -t64 * t105 + t106 * t76;
t231 = t52 * t88;
t111 = legFrame(2,3);
t65 = sin(t111);
t77 = cos(t111);
t53 = -t65 * t105 + t106 * t77;
t230 = t53 * t89;
t112 = legFrame(1,3);
t66 = sin(t112);
t78 = cos(t112);
t54 = -t66 * t105 + t106 * t78;
t229 = t54 * t90;
t23 = -t153 * t98 + t178;
t59 = t100 * t138 + t144 * t97;
t17 = t23 * t102 + t99 * t59;
t222 = t105 * t17;
t24 = -t152 * t98 + t177;
t60 = t100 * t139 + t145 * t97;
t18 = t24 * t102 + t99 * t60;
t221 = t105 * t18;
t220 = t105 * t25;
t219 = t105 * t26;
t218 = t105 * t27;
t217 = t105 * t28;
t216 = t105 * t29;
t215 = t105 * t30;
t214 = t105 * t85;
t213 = t105 * t86;
t212 = t105 * t87;
t211 = t105 * t88;
t210 = t106 * t17;
t209 = t106 * t18;
t208 = t106 * t85;
t207 = t106 * t86;
t206 = t106 * t87;
t205 = t106 * t88;
t11 = t102 * t59 - t99 * t23;
t204 = t11 * t105;
t203 = t11 * t106;
t12 = t102 * t60 - t99 * t24;
t202 = t12 * t105;
t201 = t12 * t106;
t188 = t25 * t106;
t187 = t26 * t106;
t186 = t27 * t106;
t185 = t28 * t106;
t184 = t29 * t106;
t183 = t30 * t106;
t164 = 0.1e1 / sin(pkin(9)) / sin(pkin(5));
t84 = cos(t118);
t83 = cos(t117);
t82 = cos(t116);
t81 = cos(t115);
t80 = cos(t114);
t79 = cos(t113);
t58 = t100 * t137 + t143 * t97;
t57 = t100 * t136 + t142 * t97;
t56 = t100 * t135 + t141 * t97;
t55 = t100 * t134 + t140 * t97;
t48 = t78 * t105 + t66 * t106;
t47 = t77 * t105 + t65 * t106;
t46 = t76 * t105 + t64 * t106;
t45 = t75 * t105 + t63 * t106;
t44 = t74 * t105 + t62 * t106;
t43 = t73 * t105 + t61 * t106;
t22 = -t154 * t98 + t179;
t21 = -t155 * t98 + t180;
t20 = -t156 * t98 + t181;
t19 = -t157 * t98 + t182;
t16 = t22 * t102 + t99 * t58;
t15 = t21 * t102 + t99 * t57;
t14 = t20 * t102 + t99 * t56;
t13 = t19 * t102 + t99 * t55;
t10 = t102 * t58 - t99 * t22;
t9 = t102 * t57 - t99 * t21;
t8 = t102 * t56 - t99 * t20;
t7 = t102 * t55 - t99 * t19;
t6 = (t102 * t177 + (-t98 * t165 + t195) * t100 + t97 * (t98 * t171 + t189)) * t96 + t84 * t223;
t5 = (t102 * t178 + (-t98 * t166 + t196) * t100 + t97 * (t98 * t172 + t190)) * t95 + t83 * t224;
t4 = (t102 * t179 + (-t98 * t167 + t197) * t100 + t97 * (t98 * t173 + t191)) * t94 + t82 * t225;
t3 = (t102 * t180 + (-t98 * t168 + t198) * t100 + t97 * (t98 * t174 + t192)) * t93 + t81 * t226;
t2 = (t102 * t181 + (-t98 * t169 + t199) * t100 + t97 * (t98 * t175 + t193)) * t92 + t80 * t227;
t1 = (t102 * t182 + (-t98 * t170 + t200) * t100 + t97 * (t98 * t176 + t194)) * t91 + t79 * t228;
t31 = [t54 * t96 * t164, (t72 * t229 + t48 * t84) * t164, (-t84 * t229 + t72 * t48) * t164, (((-t90 * t201 - t221) * t78 - (-t90 * t202 + t209) * t66) * t84 - t72 * ((t90 * t209 - t202) * t78 - t66 * (t90 * t221 + t201))) * t164, ((t6 * t106 - t72 * t215) * t78 - t66 * (t105 * t6 + t72 * t183)) * t164, ((t240 * t106 + t84 * t215) * t78 + t66 * (-t240 * t105 + t84 * t183)) * t164; t53 * t95 * t164, (t71 * t230 + t47 * t83) * t164, (-t83 * t230 + t71 * t47) * t164, (((-t89 * t203 - t222) * t77 - (-t89 * t204 + t210) * t65) * t83 - t71 * ((t89 * t210 - t204) * t77 - t65 * (t89 * t222 + t203))) * t164, ((t5 * t106 - t71 * t216) * t77 - t65 * (t105 * t5 + t71 * t184)) * t164, ((t239 * t106 + t83 * t216) * t77 + t65 * (-t239 * t105 + t83 * t184)) * t164; t52 * t94 * t164, (t70 * t231 + t46 * t82) * t164, (-t82 * t231 + t70 * t46) * t164, (((-t10 * t205 - t105 * t16) * t76 - (-t10 * t211 + t106 * t16) * t64) * t82 - t70 * ((-t10 * t105 + t16 * t205) * t76 - t64 * (t10 * t106 + t16 * t211))) * t164, ((t4 * t106 - t70 * t217) * t76 - t64 * (t105 * t4 + t70 * t185)) * t164, ((t238 * t106 + t82 * t217) * t76 + t64 * (-t238 * t105 + t82 * t185)) * t164; t51 * t93 * t164, (t69 * t232 + t45 * t81) * t164, (-t81 * t232 + t69 * t45) * t164, (((-t105 * t15 - t9 * t206) * t75 - (t106 * t15 - t9 * t212) * t63) * t81 - t69 * ((-t9 * t105 + t15 * t206) * t75 - t63 * (t9 * t106 + t15 * t212))) * t164, ((t3 * t106 - t69 * t218) * t75 - t63 * (t105 * t3 + t69 * t186)) * t164, ((t237 * t106 + t81 * t218) * t75 + t63 * (-t237 * t105 + t81 * t186)) * t164; t50 * t92 * t164, (t68 * t233 + t44 * t80) * t164, (-t80 * t233 + t68 * t44) * t164, (((-t105 * t14 - t8 * t207) * t74 - (t106 * t14 - t8 * t213) * t62) * t80 - t68 * ((-t8 * t105 + t14 * t207) * t74 - t62 * (t8 * t106 + t14 * t213))) * t164, ((t2 * t106 - t68 * t219) * t74 - t62 * (t105 * t2 + t68 * t187)) * t164, ((t236 * t106 + t80 * t219) * t74 + t62 * (-t236 * t105 + t80 * t187)) * t164; t49 * t91 * t164, (t67 * t234 + t43 * t79) * t164, (-t79 * t234 + t67 * t43) * t164, (((-t105 * t13 - t7 * t208) * t73 - (t106 * t13 - t7 * t214) * t61) * t79 - t67 * ((-t7 * t105 + t13 * t208) * t73 - t61 * (t7 * t106 + t13 * t214))) * t164, ((t1 * t106 - t67 * t220) * t73 - t61 * (t105 * t1 + t67 * t188)) * t164, ((t235 * t106 + t79 * t220) * t73 + t61 * (-t235 * t105 + t79 * t188)) * t164;];
Jinv  = t31;
