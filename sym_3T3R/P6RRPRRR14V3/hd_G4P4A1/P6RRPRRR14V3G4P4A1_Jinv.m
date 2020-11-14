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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2020-09-28 23:15
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P6RRPRRR14V3G4P4A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,6),zeros(1,1),zeros(6,3),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [3 6]), ...
  'P6RRPRRR14V3G4P4A1_Jinv: qJ has to be [3x6] (double)');
assert(isreal(xP) && all(size(xP) == [6 1]), ...
  'P6RRPRRR14V3G4P4A1_Jinv: xP has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'P6RRPRRR14V3G4P4A1_Jinv: pkin has to be [1x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [6 3]), ...
  'P6RRPRRR14V3G4P4A1_Jinv: legFrame has to be [6x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6RRPRRR14V3G4P4A1_Jinv: Koppelpunkt has to be [6x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-09-28 23:14:51
% EndTime: 2020-09-28 23:14:52
% DurationCPUTime: 0.95s
% Computational Cost: add. (318->159), mult. (843->330), div. (0->0), fcn. (1041->66), ass. (0->223)
t122 = cos(qJ(2,4));
t116 = sin(qJ(2,4));
t117 = sin(qJ(1,4));
t123 = cos(qJ(1,4));
t102 = legFrame(4,3);
t60 = sin(t102);
t72 = cos(t102);
t42 = -t60 * t117 + t72 * t123;
t163 = t116 * t42;
t126 = legFrame(4,2);
t84 = sin(t126);
t90 = cos(t126);
t24 = -t84 * t122 + t90 * t163;
t120 = cos(qJ(2,5));
t114 = sin(qJ(2,5));
t115 = sin(qJ(1,5));
t121 = cos(qJ(1,5));
t101 = legFrame(5,3);
t59 = sin(t101);
t71 = cos(t101);
t41 = -t59 * t115 + t71 * t121;
t164 = t114 * t41;
t125 = legFrame(5,2);
t83 = sin(t125);
t89 = cos(t125);
t23 = -t83 * t120 + t89 * t164;
t118 = cos(qJ(2,6));
t112 = sin(qJ(2,6));
t113 = sin(qJ(1,6));
t119 = cos(qJ(1,6));
t100 = legFrame(6,3);
t58 = sin(t100);
t70 = cos(t100);
t40 = -t58 * t113 + t70 * t119;
t165 = t112 * t40;
t124 = legFrame(6,2);
t82 = sin(t124);
t88 = cos(t124);
t22 = -t82 * t118 + t88 * t165;
t52 = t118 * t88;
t106 = legFrame(6,1);
t64 = sin(t106);
t228 = t52 * t64;
t53 = t120 * t89;
t107 = legFrame(5,1);
t65 = sin(t107);
t227 = t53 * t65;
t54 = t122 * t90;
t108 = legFrame(4,1);
t66 = sin(t108);
t226 = t54 * t66;
t225 = t64 * t82;
t224 = t65 * t83;
t223 = t66 * t84;
t137 = cos(qJ(1,3));
t131 = sin(qJ(1,3));
t103 = legFrame(3,3);
t73 = cos(t103);
t209 = t131 * t73;
t61 = sin(t103);
t43 = t137 * t61 + t209;
t109 = legFrame(3,1);
t79 = cos(t109);
t222 = t79 * t43;
t127 = legFrame(3,2);
t85 = sin(t127);
t221 = t79 * t85;
t139 = cos(qJ(1,2));
t133 = sin(qJ(1,2));
t104 = legFrame(2,3);
t74 = cos(t104);
t207 = t133 * t74;
t62 = sin(t104);
t44 = t139 * t62 + t207;
t110 = legFrame(2,1);
t80 = cos(t110);
t220 = t80 * t44;
t128 = legFrame(2,2);
t86 = sin(t128);
t219 = t80 * t86;
t141 = cos(qJ(1,1));
t135 = sin(qJ(1,1));
t105 = legFrame(1,3);
t75 = cos(t105);
t205 = t135 * t75;
t63 = sin(t105);
t45 = t141 * t63 + t205;
t111 = legFrame(1,1);
t81 = cos(t111);
t218 = t81 * t45;
t129 = legFrame(1,2);
t87 = sin(t129);
t217 = t81 * t87;
t176 = t73 * t137;
t130 = sin(qJ(2,3));
t210 = t130 * t85;
t136 = cos(qJ(2,3));
t91 = cos(t127);
t55 = t136 * t91;
t216 = t176 * t210 + t55;
t175 = t74 * t139;
t132 = sin(qJ(2,2));
t208 = t132 * t86;
t138 = cos(qJ(2,2));
t92 = cos(t128);
t56 = t138 * t92;
t215 = t175 * t208 + t56;
t174 = t75 * t141;
t134 = sin(qJ(2,1));
t206 = t134 * t87;
t140 = cos(qJ(2,1));
t93 = cos(t129);
t57 = t140 * t93;
t214 = t174 * t206 + t57;
t37 = t113 * t70 + t119 * t58;
t213 = t112 * t37;
t38 = t115 * t71 + t121 * t59;
t212 = t114 * t38;
t39 = t117 * t72 + t123 * t60;
t211 = t116 * t39;
t150 = koppelP(1,3);
t143 = xP(5);
t98 = cos(t143);
t204 = t150 * t98;
t151 = koppelP(6,2);
t144 = xP(4);
t96 = sin(t144);
t203 = t151 * t96;
t99 = cos(t144);
t202 = t151 * t99;
t152 = koppelP(5,2);
t201 = t152 * t96;
t200 = t152 * t99;
t153 = koppelP(4,2);
t199 = t153 * t96;
t198 = t153 * t99;
t154 = koppelP(3,2);
t197 = t154 * t96;
t196 = t154 * t99;
t155 = koppelP(2,2);
t195 = t155 * t96;
t194 = t155 * t99;
t156 = koppelP(1,2);
t193 = t156 * t96;
t192 = t156 * t99;
t157 = koppelP(6,1);
t191 = t157 * t96;
t190 = t157 * t99;
t158 = koppelP(5,1);
t189 = t158 * t96;
t188 = t158 * t99;
t159 = koppelP(4,1);
t187 = t159 * t96;
t186 = t159 * t99;
t160 = koppelP(3,1);
t185 = t160 * t96;
t184 = t160 * t99;
t161 = koppelP(2,1);
t183 = t161 * t96;
t182 = t161 * t99;
t162 = koppelP(1,1);
t181 = t162 * t96;
t180 = t162 * t99;
t179 = t61 * t131;
t178 = t62 * t133;
t177 = t63 * t135;
t145 = koppelP(6,3);
t170 = t98 * t145;
t146 = koppelP(5,3);
t169 = t98 * t146;
t147 = koppelP(4,3);
t168 = t98 * t147;
t148 = koppelP(3,3);
t167 = t98 * t148;
t149 = koppelP(2,3);
t166 = t98 * t149;
t142 = xP(6);
t97 = cos(t142);
t95 = sin(t143);
t94 = sin(t142);
t78 = cos(t108);
t77 = cos(t107);
t76 = cos(t106);
t69 = sin(t111);
t68 = sin(t110);
t67 = sin(t109);
t48 = t174 - t177;
t47 = t175 - t178;
t46 = t176 - t179;
t36 = (-t156 * t94 + t162 * t97) * t98 + t95 * t150;
t35 = (-t155 * t94 + t161 * t97) * t98 + t95 * t149;
t34 = (-t154 * t94 + t160 * t97) * t98 + t95 * t148;
t33 = (-t153 * t94 + t159 * t97) * t98 + t95 * t147;
t32 = (-t152 * t94 + t158 * t97) * t98 + t95 * t146;
t31 = (-t151 * t94 + t157 * t97) * t98 + t95 * t145;
t30 = t134 * t48 * t93 - t87 * t140;
t29 = t132 * t47 * t92 - t86 * t138;
t28 = t130 * t46 * t91 - t85 * t136;
t27 = t84 * t163 + t54;
t26 = t83 * t164 + t53;
t25 = t82 * t165 + t52;
t21 = (-t95 * t192 - t181) * t94 + (t95 * t180 - t193) * t97 - t99 * t204;
t20 = (-t95 * t193 + t180) * t94 + (t95 * t181 + t192) * t97 - t96 * t204;
t19 = (t95 * t182 - t195) * t97 + (-t95 * t194 - t183) * t94 - t99 * t166;
t18 = (t95 * t184 - t197) * t97 + (-t95 * t196 - t185) * t94 - t99 * t167;
t17 = (t95 * t186 - t199) * t97 + (-t95 * t198 - t187) * t94 - t99 * t168;
t16 = (t95 * t188 - t201) * t97 + (-t95 * t200 - t189) * t94 - t99 * t169;
t15 = (t95 * t190 - t203) * t97 + (-t95 * t202 - t191) * t94 - t99 * t170;
t14 = -t96 * t166 + (t95 * t183 + t194) * t97 + t94 * (-t95 * t195 + t182);
t13 = -t96 * t167 + (t95 * t185 + t196) * t97 + t94 * (-t95 * t197 + t184);
t12 = -t96 * t168 + (t95 * t187 + t198) * t97 + t94 * (-t95 * t199 + t186);
t11 = -t96 * t169 + (t95 * t189 + t200) * t97 + t94 * (-t95 * t201 + t188);
t10 = -t96 * t170 + (t95 * t191 + t202) * t97 + t94 * (-t95 * t203 + t190);
t9 = (-t177 * t206 + t214) * t69 + t134 * t218;
t8 = (-t178 * t208 + t215) * t68 + t132 * t220;
t7 = (-t179 * t210 + t216) * t67 + t130 * t222;
t6 = -t66 * t211 + t27 * t78;
t5 = -t65 * t212 + t26 * t77;
t4 = -t64 * t213 + t25 * t76;
t3 = t81 * t214 + (-(t135 * t217 + t141 * t69) * t63 - t69 * t205) * t134;
t2 = t80 * t215 + (-(t133 * t219 + t68 * t139) * t62 - t68 * t207) * t132;
t1 = t79 * t216 + (-(t131 * t221 + t67 * t137) * t61 - t67 * t209) * t130;
t49 = [t30, (t87 * t48 * t69 + t218) * t134 + t69 * t57, (-t48 * t217 + t69 * t45) * t134 - t81 * t57, -t20 * t3 + t9 * t21, -t21 * t30 + t3 * t36, -t30 * t20 + t9 * t36; t29, (t86 * t47 * t68 + t220) * t132 + t68 * t56, (-t47 * t219 + t68 * t44) * t132 - t80 * t56, -t14 * t2 + t8 * t19, -t19 * t29 + t2 * t35, -t29 * t14 + t8 * t35; t28, (t85 * t46 * t67 + t222) * t130 + t67 * t55, (-t46 * t221 + t67 * t43) * t130 - t79 * t55, -t13 * t1 + t7 * t18, t1 * t34 - t18 * t28, -t28 * t13 + t7 * t34; t24, (t42 * t223 + t78 * t39) * t116 + t226, (-t84 * t42 * t78 + t66 * t39) * t116 - t78 * t54, (t226 + ((t117 * t78 + t123 * t223) * t72 + (-t117 * t223 + t123 * t78) * t60) * t116) * t17 - t12 * t6, -t24 * t17 + t33 * t6, -t24 * t12 + t33 * (t78 * t211 + t27 * t66); t23, (t41 * t224 + t77 * t38) * t114 + t227, (-t83 * t41 * t77 + t65 * t38) * t114 - t77 * t53, (t227 + ((t115 * t77 + t121 * t224) * t71 + (-t115 * t224 + t121 * t77) * t59) * t114) * t16 - t11 * t5, -t23 * t16 + t32 * t5, -t23 * t11 + t32 * (t77 * t212 + t26 * t65); t22, (t40 * t225 + t76 * t37) * t112 + t228, (-t82 * t40 * t76 + t64 * t37) * t112 - t76 * t52, (t228 + ((t113 * t76 + t119 * t225) * t70 + (-t113 * t225 + t119 * t76) * t58) * t112) * t15 - t10 * t4, -t22 * t15 + t31 * t4, -t22 * t10 + t31 * (t76 * t213 + t25 * t64);];
Jinv  = t49;
