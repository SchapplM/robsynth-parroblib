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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,alpha3,alpha4,d1,d4,theta3]';
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
% Datum: 2022-11-03 12:08
% Revision: 5de314b2d97380370c92dd342c670b987640c890 (2022-02-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P6RRPRRR14V4G7P4A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,6),zeros(7,1),zeros(6,3),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [3 6]), ...
  'P6RRPRRR14V4G7P4A1_Jinv: qJ has to be [3x6] (double)');
assert(isreal(xP) && all(size(xP) == [6 1]), ...
  'P6RRPRRR14V4G7P4A1_Jinv: xP has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P6RRPRRR14V4G7P4A1_Jinv: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [6 3]), ...
  'P6RRPRRR14V4G7P4A1_Jinv: legFrame has to be [6x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6RRPRRR14V4G7P4A1_Jinv: Koppelpunkt has to be [6x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-03 12:08:35
% EndTime: 2022-11-03 12:08:37
% DurationCPUTime: 2.15s
% Computational Cost: add. (984->273), mult. (2244->548), div. (54->6), fcn. (2448->48), ass. (0->262)
t157 = xP(5);
t115 = sin(t157);
t118 = cos(t157);
t164 = koppelP(1,3);
t156 = xP(6);
t114 = sin(t156);
t117 = cos(t156);
t170 = koppelP(1,2);
t176 = koppelP(1,1);
t189 = t114 * t170 - t117 * t176;
t305 = -t115 * t164 + t118 * t189;
t163 = koppelP(2,3);
t169 = koppelP(2,2);
t175 = koppelP(2,1);
t190 = t114 * t169 - t117 * t175;
t304 = -t115 * t163 + t118 * t190;
t162 = koppelP(3,3);
t168 = koppelP(3,2);
t174 = koppelP(3,1);
t191 = t114 * t168 - t117 * t174;
t303 = -t115 * t162 + t118 * t191;
t161 = koppelP(4,3);
t167 = koppelP(4,2);
t173 = koppelP(4,1);
t192 = t114 * t167 - t117 * t173;
t302 = -t115 * t161 + t118 * t192;
t160 = koppelP(5,3);
t166 = koppelP(5,2);
t172 = koppelP(5,1);
t193 = t114 * t166 - t117 * t172;
t301 = -t115 * t160 + t118 * t193;
t159 = koppelP(6,3);
t165 = koppelP(6,2);
t171 = koppelP(6,1);
t194 = t114 * t165 - t117 * t171;
t300 = -t115 * t159 + t118 * t194;
t126 = sin(qJ(2,6));
t132 = cos(qJ(2,6));
t122 = sin(pkin(3));
t125 = cos(pkin(3));
t120 = sin(pkin(7));
t123 = cos(pkin(7));
t121 = sin(pkin(4));
t281 = pkin(6) * t121;
t88 = -t120 * pkin(2) + t123 * t281;
t262 = t88 * t125;
t124 = cos(pkin(4));
t101 = t124 * pkin(6);
t95 = t101 + qJ(3,6);
t288 = t122 * t95 + t262;
t81 = t123 * pkin(2) + t120 * t281 + pkin(1);
t299 = t288 * t126 + t132 * t81;
t128 = sin(qJ(2,5));
t134 = cos(qJ(2,5));
t96 = t101 + qJ(3,5);
t289 = t122 * t96 + t262;
t298 = t289 * t128 + t134 * t81;
t130 = sin(qJ(2,4));
t136 = cos(qJ(2,4));
t97 = t101 + qJ(3,4);
t290 = t122 * t97 + t262;
t297 = t290 * t130 + t136 * t81;
t144 = sin(qJ(2,3));
t150 = cos(qJ(2,3));
t98 = t101 + qJ(3,3);
t291 = t122 * t98 + t262;
t296 = t291 * t144 + t150 * t81;
t146 = sin(qJ(2,2));
t152 = cos(qJ(2,2));
t99 = t101 + qJ(3,2);
t292 = t122 * t99 + t262;
t295 = t292 * t146 + t152 * t81;
t148 = sin(qJ(2,1));
t154 = cos(qJ(2,1));
t100 = t101 + qJ(3,1);
t293 = t122 * t100 + t262;
t294 = t293 * t148 + t154 * t81;
t127 = sin(qJ(1,6));
t133 = cos(qJ(1,6));
t280 = t122 * t88;
t46 = t95 * t125 - t280;
t287 = t127 * t46 + t299 * t133;
t129 = sin(qJ(1,5));
t135 = cos(qJ(1,5));
t47 = t96 * t125 - t280;
t286 = t129 * t47 + t298 * t135;
t131 = sin(qJ(1,4));
t137 = cos(qJ(1,4));
t48 = t97 * t125 - t280;
t285 = t131 * t48 + t297 * t137;
t145 = sin(qJ(1,3));
t151 = cos(qJ(1,3));
t52 = t98 * t125 - t280;
t284 = t145 * t52 + t296 * t151;
t147 = sin(qJ(1,2));
t153 = cos(qJ(1,2));
t53 = t99 * t125 - t280;
t283 = t147 * t53 + t295 * t153;
t149 = sin(qJ(1,1));
t155 = cos(qJ(1,1));
t54 = t100 * t125 - t280;
t282 = t149 * t54 + t294 * t155;
t158 = xP(4);
t116 = sin(t158);
t255 = t116 * t165;
t254 = t116 * t166;
t253 = t116 * t167;
t252 = t116 * t168;
t251 = t116 * t169;
t250 = t116 * t170;
t249 = t116 * t171;
t248 = t116 * t172;
t247 = t116 * t173;
t246 = t116 * t174;
t245 = t116 * t175;
t244 = t116 * t176;
t243 = t118 * t159;
t242 = t118 * t160;
t241 = t118 * t161;
t240 = t118 * t162;
t239 = t118 * t163;
t238 = t118 * t164;
t119 = cos(t158);
t237 = t119 * t165;
t236 = t119 * t166;
t235 = t119 * t167;
t234 = t119 * t168;
t233 = t119 * t169;
t232 = t119 * t170;
t231 = t119 * t171;
t230 = t119 * t172;
t229 = t119 * t173;
t228 = t119 * t174;
t227 = t119 * t175;
t226 = t119 * t176;
t224 = t122 * t127;
t223 = t122 * t129;
t222 = t122 * t131;
t221 = t122 * t145;
t220 = t122 * t147;
t219 = t122 * t149;
t218 = t125 * t127;
t217 = t125 * t129;
t216 = t125 * t131;
t215 = t125 * t145;
t214 = t125 * t147;
t213 = t125 * t149;
t212 = t126 * t133;
t211 = t128 * t135;
t210 = t130 * t137;
t209 = t132 * t133;
t208 = t134 * t135;
t207 = t136 * t137;
t206 = t144 * t151;
t205 = t146 * t153;
t204 = t148 * t155;
t203 = t150 * t151;
t202 = t152 * t153;
t201 = t154 * t155;
t138 = legFrame(6,2);
t102 = sin(t138);
t108 = cos(t138);
t55 = t102 * t212 + t108 * t132;
t188 = t102 * t218 + t122 * t55;
t139 = legFrame(5,2);
t103 = sin(t139);
t109 = cos(t139);
t56 = t103 * t211 + t109 * t134;
t187 = t103 * t217 + t122 * t56;
t140 = legFrame(4,2);
t104 = sin(t140);
t110 = cos(t140);
t57 = t104 * t210 + t110 * t136;
t186 = t104 * t216 + t122 * t57;
t141 = legFrame(3,2);
t105 = sin(t141);
t111 = cos(t141);
t67 = t105 * t206 + t111 * t150;
t185 = t105 * t215 + t122 * t67;
t142 = legFrame(2,2);
t106 = sin(t142);
t112 = cos(t142);
t68 = t106 * t205 + t112 * t152;
t184 = t106 * t214 + t122 * t68;
t143 = legFrame(1,2);
t107 = sin(t143);
t113 = cos(t143);
t69 = t107 * t204 + t113 * t154;
t183 = t107 * t213 + t122 * t69;
t58 = -t132 * t102 + t108 * t212;
t182 = t108 * t218 + t58 * t122;
t59 = -t134 * t103 + t109 * t211;
t181 = t109 * t217 + t59 * t122;
t60 = -t136 * t104 + t110 * t210;
t180 = t110 * t216 + t60 * t122;
t71 = -t150 * t105 + t111 * t206;
t179 = t111 * t215 + t71 * t122;
t73 = -t152 * t106 + t112 * t205;
t178 = t112 * t214 + t73 * t122;
t74 = -t154 * t107 + t113 * t204;
t177 = t113 * t213 + t74 * t122;
t94 = 0.1e1 / t100;
t93 = 0.1e1 / t99;
t92 = 0.1e1 / t98;
t91 = 0.1e1 / t97;
t90 = 0.1e1 / t96;
t89 = 0.1e1 / t95;
t78 = t106 * t202 - t112 * t146;
t77 = t105 * t203 - t144 * t111;
t76 = t144 * t105 + t111 * t203;
t75 = t148 * t107 + t113 * t201;
t72 = t146 * t106 + t112 * t202;
t70 = t107 * t201 - t148 * t113;
t66 = t104 * t207 - t130 * t110;
t65 = t130 * t104 + t110 * t207;
t64 = t103 * t208 - t128 * t109;
t63 = t128 * t103 + t109 * t208;
t62 = t102 * t209 - t126 * t108;
t61 = t126 * t102 + t108 * t209;
t42 = -t113 * t219 + t74 * t125;
t41 = -t112 * t220 + t73 * t125;
t40 = -t111 * t221 + t71 * t125;
t39 = -t107 * t219 + t69 * t125;
t38 = -t106 * t220 + t68 * t125;
t37 = -t105 * t221 + t67 * t125;
t36 = -t110 * t222 + t60 * t125;
t35 = -t109 * t223 + t59 * t125;
t34 = -t108 * t224 + t58 * t125;
t33 = -t104 * t222 + t57 * t125;
t32 = -t103 * t223 + t56 * t125;
t31 = -t102 * t224 + t55 * t125;
t30 = t81 * t148 - t293 * t154;
t29 = t81 * t146 - t292 * t152;
t28 = t81 * t144 - t291 * t150;
t27 = t81 * t130 - t290 * t136;
t26 = t81 * t128 - t289 * t134;
t25 = t81 * t126 - t288 * t132;
t24 = (t115 * t226 - t250) * t117 + (-t115 * t232 - t244) * t114 - t119 * t238;
t23 = (t115 * t227 - t251) * t117 + (-t115 * t233 - t245) * t114 - t119 * t239;
t22 = (t115 * t228 - t252) * t117 + (-t115 * t234 - t246) * t114 - t119 * t240;
t21 = (t115 * t229 - t253) * t117 + (-t115 * t235 - t247) * t114 - t119 * t241;
t20 = (t115 * t230 - t254) * t117 + (-t115 * t236 - t248) * t114 - t119 * t242;
t19 = (t115 * t231 - t255) * t117 + (-t115 * t237 - t249) * t114 - t119 * t243;
t18 = -t294 * t149 + t54 * t155;
t17 = -t295 * t147 + t53 * t153;
t16 = -t296 * t145 + t52 * t151;
t15 = -t297 * t131 + t48 * t137;
t14 = -t298 * t129 + t47 * t135;
t13 = -t299 * t127 + t46 * t133;
t12 = pkin(1) * t78 + t184 * qJ(3,2) + (-t38 * t120 + t78 * t123) * pkin(2) + (t184 * t124 + (t120 * t78 + t38 * t123) * t121) * pkin(6);
t11 = pkin(1) * t77 + t185 * qJ(3,3) + (-t37 * t120 + t77 * t123) * pkin(2) + (t185 * t124 + (t120 * t77 + t37 * t123) * t121) * pkin(6);
t10 = pkin(1) * t75 + t177 * qJ(3,1) + (-t42 * t120 + t75 * t123) * pkin(2) + (t177 * t124 + (t120 * t75 + t42 * t123) * t121) * pkin(6);
t9 = pkin(1) * t72 + t178 * qJ(3,2) + (-t41 * t120 + t72 * t123) * pkin(2) + (t178 * t124 + (t120 * t72 + t41 * t123) * t121) * pkin(6);
t8 = pkin(1) * t76 + t179 * qJ(3,3) + (-t40 * t120 + t76 * t123) * pkin(2) + (t179 * t124 + (t120 * t76 + t40 * t123) * t121) * pkin(6);
t7 = pkin(1) * t70 + t183 * qJ(3,1) + (-t39 * t120 + t70 * t123) * pkin(2) + (t183 * t124 + (t120 * t70 + t39 * t123) * t121) * pkin(6);
t6 = pkin(1) * t65 + t180 * qJ(3,4) + (-t36 * t120 + t65 * t123) * pkin(2) + (t180 * t124 + (t120 * t65 + t36 * t123) * t121) * pkin(6);
t5 = pkin(1) * t63 + t181 * qJ(3,5) + (-t35 * t120 + t63 * t123) * pkin(2) + (t181 * t124 + (t120 * t63 + t35 * t123) * t121) * pkin(6);
t4 = pkin(1) * t61 + t182 * qJ(3,6) + (-t34 * t120 + t61 * t123) * pkin(2) + (t182 * t124 + (t120 * t61 + t34 * t123) * t121) * pkin(6);
t3 = pkin(1) * t66 + t186 * qJ(3,4) + (-t33 * t120 + t66 * t123) * pkin(2) + (t186 * t124 + (t120 * t66 + t33 * t123) * t121) * pkin(6);
t2 = pkin(1) * t64 + t187 * qJ(3,5) + (-t32 * t120 + t64 * t123) * pkin(2) + (t187 * t124 + (t120 * t64 + t32 * t123) * t121) * pkin(6);
t1 = pkin(1) * t62 + t188 * qJ(3,6) + (-t31 * t120 + t62 * t123) * pkin(2) + (t188 * t124 + (t120 * t62 + t31 * t123) * t121) * pkin(6);
t43 = [(t30 * t107 + t282 * t113) * t94, (-t282 * t107 + t113 * t30) * t94, t18 * t94, (-t24 * t7 + t18 * ((-t115 * t189 - t238) * t116 + t119 * (t114 * t176 + t117 * t170))) * t94, (-t24 * t10 + t18 * t305) * t94, (-(-t116 * t238 + (t115 * t244 + t232) * t117 + t114 * (-t115 * t250 + t226)) * t10 + t305 * t7) * t94; (t29 * t106 + t283 * t112) * t93, (-t283 * t106 + t112 * t29) * t93, t17 * t93, (-t23 * t12 + t17 * ((-t115 * t190 - t239) * t116 + t119 * (t114 * t175 + t117 * t169))) * t93, (t17 * t304 - t23 * t9) * t93, (-(-t116 * t239 + (t115 * t245 + t233) * t117 + t114 * (-t115 * t251 + t227)) * t9 + t304 * t12) * t93; (t28 * t105 + t284 * t111) * t92, (-t284 * t105 + t111 * t28) * t92, t16 * t92, (-t22 * t11 + t16 * ((-t115 * t191 - t240) * t116 + t119 * (t114 * t174 + t117 * t168))) * t92, (t16 * t303 - t22 * t8) * t92, (-(-t116 * t240 + (t115 * t246 + t234) * t117 + t114 * (-t115 * t252 + t228)) * t8 + t303 * t11) * t92; (t27 * t104 + t285 * t110) * t91, (-t285 * t104 + t110 * t27) * t91, t15 * t91, (-t21 * t3 + t15 * ((-t115 * t192 - t241) * t116 + t119 * (t114 * t173 + t117 * t167))) * t91, (t15 * t302 - t21 * t6) * t91, (-t6 * (-t116 * t241 + (t115 * t247 + t235) * t117 + t114 * (-t115 * t253 + t229)) + t302 * t3) * t91; (t26 * t103 + t286 * t109) * t90, (-t286 * t103 + t109 * t26) * t90, t14 * t90, (-t20 * t2 + t14 * ((-t115 * t193 - t242) * t116 + t119 * (t114 * t172 + t117 * t166))) * t90, (t14 * t301 - t20 * t5) * t90, (-t5 * (-t116 * t242 + (t115 * t248 + t236) * t117 + t114 * (-t115 * t254 + t230)) + t301 * t2) * t90; (t25 * t102 + t287 * t108) * t89, (-t287 * t102 + t108 * t25) * t89, t13 * t89, (-t19 * t1 + t13 * ((-t115 * t194 - t243) * t116 + t119 * (t114 * t171 + t117 * t165))) * t89, (t13 * t300 - t19 * t4) * t89, (-t4 * (-t116 * t243 + (t115 * t249 + t237) * t117 + t114 * (-t115 * t255 + t231)) + t300 * t1) * t89;];
Jinv  = t43;
