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
% Datum: 2019-05-16 19:48
% Revision: 36f6366a01c4a552c0708fcd8ed3e0fb9da693e2 (2019-05-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P6PPPRRR1V2G1P1A2_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,6),zeros(10,1),zeros(6,3),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [3 6]), ...
  'P6PPPRRR1V2G1P1A2_Jinv: qJ has to be [3x6] (double)');
assert(isreal(xP) && all(size(xP) == [6 1]), ...
  'P6PPPRRR1V2G1P1A2_Jinv: xP has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'P6PPPRRR1V2G1P1A2_Jinv: pkin has to be [10x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [6 3]), ...
  'P6PPPRRR1V2G1P1A2_Jinv: legFrame has to be [6x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6PPPRRR1V2G1P1A2_Jinv: Koppelpunkt has to be [6x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-16 19:47:58
% EndTime: 2019-05-16 19:48:02
% DurationCPUTime: 3.94s
% Computational Cost: add. (1368->386), mult. (3414->869), div. (108->3), fcn. (3798->50), ass. (0->308)
t169 = legFrame(6,2);
t130 = sin(t169);
t151 = sin(pkin(9));
t152 = sin(pkin(8));
t155 = cos(pkin(8));
t232 = cos(pkin(4)) * cos(pkin(9));
t79 = t155 * t151 + t152 * t232;
t324 = t130 * t79;
t80 = t151 * t152 - t155 * t232;
t323 = t130 * t80;
t170 = legFrame(5,2);
t131 = sin(t170);
t322 = t131 * t79;
t321 = t131 * t80;
t171 = legFrame(4,2);
t132 = sin(t171);
t320 = t132 * t79;
t319 = t132 * t80;
t172 = legFrame(3,2);
t133 = sin(t172);
t318 = t133 * t79;
t317 = t133 * t80;
t173 = legFrame(2,2);
t134 = sin(t173);
t316 = t134 * t79;
t315 = t134 * t80;
t174 = legFrame(1,2);
t135 = sin(t174);
t314 = t135 * t79;
t313 = t135 * t80;
t176 = xP(5);
t143 = sin(t176);
t146 = cos(t176);
t178 = koppelP(6,3);
t175 = xP(6);
t142 = sin(t175);
t145 = cos(t175);
t184 = koppelP(6,2);
t190 = koppelP(6,1);
t219 = t142 * t184 - t145 * t190;
t73 = t143 * t178 - t219 * t146;
t312 = t152 * t73;
t179 = koppelP(5,3);
t185 = koppelP(5,2);
t191 = koppelP(5,1);
t218 = t142 * t185 - t145 * t191;
t74 = t143 * t179 - t218 * t146;
t311 = t152 * t74;
t180 = koppelP(4,3);
t186 = koppelP(4,2);
t192 = koppelP(4,1);
t217 = t142 * t186 - t145 * t192;
t75 = t143 * t180 - t217 * t146;
t310 = t152 * t75;
t181 = koppelP(3,3);
t187 = koppelP(3,2);
t193 = koppelP(3,1);
t216 = t142 * t187 - t145 * t193;
t76 = t143 * t181 - t216 * t146;
t309 = t152 * t76;
t182 = koppelP(2,3);
t188 = koppelP(2,2);
t194 = koppelP(2,1);
t215 = t142 * t188 - t145 * t194;
t77 = t143 * t182 - t215 * t146;
t308 = t152 * t77;
t183 = koppelP(1,3);
t189 = koppelP(1,2);
t195 = koppelP(1,1);
t214 = t142 * t189 - t145 * t195;
t78 = t143 * t183 - t214 * t146;
t307 = t152 * t78;
t306 = t155 * t73;
t305 = t155 * t74;
t304 = t155 * t75;
t303 = t155 * t76;
t302 = t155 * t77;
t301 = t155 * t78;
t163 = legFrame(6,1);
t112 = sin(t163);
t300 = t112 * t130;
t164 = legFrame(5,1);
t113 = sin(t164);
t299 = t113 * t131;
t165 = legFrame(4,1);
t114 = sin(t165);
t298 = t114 * t132;
t166 = legFrame(3,1);
t115 = sin(t166);
t297 = t115 * t133;
t167 = legFrame(2,1);
t116 = sin(t167);
t296 = t116 * t134;
t168 = legFrame(1,1);
t117 = sin(t168);
t295 = t117 * t135;
t124 = cos(t163);
t294 = t124 * t130;
t125 = cos(t164);
t293 = t125 * t131;
t126 = cos(t165);
t292 = t126 * t132;
t127 = cos(t166);
t291 = t127 * t133;
t128 = cos(t167);
t290 = t128 * t134;
t129 = cos(t168);
t289 = t129 * t135;
t288 = t130 * t151;
t287 = t130 * t152;
t286 = t130 * t155;
t285 = t131 * t151;
t284 = t131 * t152;
t283 = t131 * t155;
t282 = t132 * t151;
t281 = t132 * t152;
t280 = t132 * t155;
t279 = t133 * t151;
t278 = t133 * t152;
t277 = t133 * t155;
t276 = t134 * t151;
t275 = t134 * t152;
t274 = t134 * t155;
t273 = t135 * t151;
t272 = t135 * t152;
t271 = t135 * t155;
t136 = cos(t169);
t270 = t136 * t151;
t137 = cos(t170);
t269 = t137 * t151;
t138 = cos(t171);
t268 = t138 * t151;
t139 = cos(t172);
t267 = t139 * t151;
t140 = cos(t173);
t266 = t140 * t151;
t141 = cos(t174);
t265 = t141 * t151;
t177 = xP(4);
t144 = sin(t177);
t264 = t144 * t184;
t263 = t144 * t185;
t262 = t144 * t186;
t261 = t144 * t187;
t260 = t144 * t188;
t259 = t144 * t189;
t258 = t144 * t190;
t257 = t144 * t191;
t256 = t144 * t192;
t255 = t144 * t193;
t254 = t144 * t194;
t253 = t144 * t195;
t252 = t146 * t183;
t147 = cos(t177);
t251 = t147 * t184;
t250 = t147 * t185;
t249 = t147 * t186;
t248 = t147 * t187;
t247 = t147 * t188;
t246 = t147 * t189;
t245 = t147 * t190;
t244 = t147 * t191;
t243 = t147 * t192;
t242 = t147 * t193;
t241 = t147 * t194;
t240 = t147 * t195;
t154 = sin(pkin(4));
t239 = cos(pkin(5)) * t154;
t238 = t178 * t146;
t237 = t179 * t146;
t236 = t180 * t146;
t235 = t181 * t146;
t234 = t182 * t146;
t153 = sin(pkin(5));
t233 = 0.1e1 / t151 / t153 / t154;
t157 = legFrame(6,3);
t106 = sin(t157);
t118 = cos(t157);
t231 = t106 * t80 - t118 * t79;
t230 = t106 * t79 + t118 * t80;
t158 = legFrame(5,3);
t107 = sin(t158);
t119 = cos(t158);
t229 = t107 * t80 - t119 * t79;
t228 = t107 * t79 + t119 * t80;
t159 = legFrame(4,3);
t108 = sin(t159);
t120 = cos(t159);
t227 = t108 * t80 - t120 * t79;
t226 = t108 * t79 + t120 * t80;
t160 = legFrame(3,3);
t109 = sin(t160);
t121 = cos(t160);
t225 = t80 * t109 - t79 * t121;
t224 = t109 * t79 + t121 * t80;
t161 = legFrame(2,3);
t110 = sin(t161);
t122 = cos(t161);
t223 = t80 * t110 - t79 * t122;
t222 = t110 * t79 + t122 * t80;
t162 = legFrame(1,3);
t111 = sin(t162);
t123 = cos(t162);
t221 = t80 * t111 - t79 * t123;
t220 = t111 * t79 + t123 * t80;
t213 = t136 * t232;
t212 = t137 * t232;
t211 = t138 * t232;
t210 = t139 * t232;
t209 = t140 * t232;
t208 = t141 * t232;
t207 = t112 * t232;
t206 = t113 * t232;
t205 = t114 * t232;
t204 = t115 * t232;
t203 = t116 * t232;
t202 = t117 * t232;
t201 = t124 * t232;
t200 = t125 * t232;
t199 = t126 * t232;
t198 = t127 * t232;
t197 = t128 * t232;
t196 = t129 * t232;
t98 = t142 * t195 + t189 * t145;
t97 = t142 * t194 + t188 * t145;
t96 = t142 * t193 + t187 * t145;
t95 = t142 * t192 + t186 * t145;
t94 = t142 * t191 + t185 * t145;
t93 = t142 * t190 + t184 * t145;
t92 = t111 * t155 + t152 * t123;
t91 = t110 * t155 + t152 * t122;
t90 = t109 * t155 + t152 * t121;
t89 = t108 * t155 + t152 * t120;
t88 = t107 * t155 + t152 * t119;
t87 = t106 * t155 + t152 * t118;
t86 = t111 * t152 - t155 * t123;
t85 = t110 * t152 - t155 * t122;
t84 = t109 * t152 - t155 * t121;
t83 = t152 * t108 - t155 * t120;
t82 = t107 * t152 - t155 * t119;
t81 = t106 * t152 - t155 * t118;
t72 = t214 * t143 + t252;
t71 = t215 * t143 + t234;
t70 = t216 * t143 + t235;
t69 = t217 * t143 + t236;
t68 = t218 * t143 + t237;
t67 = t219 * t143 + t238;
t60 = t144 * t98 + t147 * t72;
t59 = t144 * t97 + t147 * t71;
t58 = t144 * t96 + t147 * t70;
t57 = t144 * t95 + t147 * t69;
t56 = t144 * t94 + t147 * t68;
t55 = t144 * t93 + t147 * t67;
t54 = -t144 * t72 + t147 * t98;
t53 = -t144 * t71 + t147 * t97;
t52 = -t144 * t70 + t147 * t96;
t51 = -t144 * t69 + t147 * t95;
t50 = -t144 * t68 + t147 * t94;
t49 = -t144 * t67 + t147 * t93;
t48 = t147 * t252 + (-t143 * t240 + t259) * t145 + t142 * (t143 * t246 + t253);
t47 = t147 * t234 + (-t143 * t241 + t260) * t145 + t142 * (t143 * t247 + t254);
t46 = t147 * t235 + (-t143 * t242 + t261) * t145 + t142 * (t143 * t248 + t255);
t45 = t147 * t236 + (-t143 * t243 + t262) * t145 + t142 * (t143 * t249 + t256);
t44 = t147 * t237 + (-t143 * t244 + t263) * t145 + t142 * (t143 * t250 + t257);
t43 = t147 * t238 + (-t143 * t245 + t264) * t145 + t142 * (t143 * t251 + t258);
t42 = -t144 * t252 + (t143 * t253 + t246) * t145 + t142 * (-t143 * t259 + t240);
t41 = -t144 * t234 + (t143 * t254 + t247) * t145 + t142 * (-t143 * t260 + t241);
t40 = -t144 * t235 + (t143 * t255 + t248) * t145 + t142 * (-t143 * t261 + t242);
t39 = -t144 * t236 + (t143 * t256 + t249) * t145 + t142 * (-t143 * t262 + t243);
t38 = -t144 * t237 + (t143 * t257 + t250) * t145 + t142 * (-t143 * t263 + t244);
t37 = -t144 * t238 + (t143 * t258 + t251) * t145 + t142 * (-t143 * t264 + t245);
t36 = t48 * t141 + t78 * t289;
t35 = t47 * t140 + t77 * t290;
t34 = t46 * t139 + t76 * t291;
t33 = t45 * t138 + t75 * t292;
t32 = t44 * t137 + t74 * t293;
t31 = t43 * t136 + t73 * t294;
t30 = t42 * t141 - t78 * t295;
t29 = t41 * t140 - t77 * t296;
t28 = t40 * t139 - t76 * t297;
t27 = t39 * t138 - t75 * t298;
t26 = t38 * t137 - t74 * t299;
t25 = t37 * t136 - t73 * t300;
t24 = t48 * t265 + t78 * (t129 * t273 + t202);
t23 = t47 * t266 + t77 * (t128 * t276 + t203);
t22 = t46 * t267 + t76 * (t127 * t279 + t204);
t21 = t45 * t268 + t75 * (t126 * t282 + t205);
t20 = t44 * t269 + t74 * (t125 * t285 + t206);
t19 = t43 * t270 + t73 * (t124 * t288 + t207);
t18 = t42 * t265 + t78 * (-t117 * t273 + t196);
t17 = t41 * t266 + t77 * (-t116 * t276 + t197);
t16 = t40 * t267 + t76 * (-t115 * t279 + t198);
t15 = t39 * t268 + t75 * (-t114 * t282 + t199);
t14 = t38 * t269 + t74 * (-t113 * t285 + t200);
t13 = t37 * t270 + t73 * (-t112 * t288 + t201);
t12 = t48 * t208 + t78 * (-t117 * t151 + t135 * t196);
t11 = t47 * t209 + t77 * (-t116 * t151 + t134 * t197);
t10 = t46 * t210 + t76 * (-t115 * t151 + t133 * t198);
t9 = t45 * t211 + t75 * (-t114 * t151 + t132 * t199);
t8 = t44 * t212 + t74 * (-t113 * t151 + t131 * t200);
t7 = t43 * t213 + t73 * (-t112 * t151 + t130 * t201);
t6 = t42 * t208 - t78 * (t129 * t151 + t135 * t202);
t5 = t41 * t209 - t77 * (t128 * t151 + t134 * t203);
t4 = t40 * t210 - t76 * (t127 * t151 + t133 * t204);
t3 = t39 * t211 - t75 * (t126 * t151 + t132 * t205);
t2 = t38 * t212 - t74 * (t125 * t151 + t131 * t206);
t1 = t37 * t213 - t73 * (t124 * t151 + t130 * t207);
t61 = [t141 * (t220 * t153 + t86 * t239) * t233, ((t221 * t129 + t220 * t295) * t153 - (t129 * t92 - t86 * t295) * t239) * t233, ((t221 * t117 - t220 * t289) * t153 + (-t117 * t92 - t86 * t289) * t239) * t233, ((((-t54 * t313 + t60 * t79) * t123 + t111 * (-t54 * t314 - t60 * t80)) * t129 + t117 * ((-t60 * t313 - t54 * t79) * t123 - (t60 * t314 - t54 * t80) * t111)) * t153 + (((t152 * t60 + t54 * t271) * t123 + t111 * (t155 * t60 - t54 * t272)) * t129 + t117 * ((-t152 * t54 + t60 * t271) * t123 - t111 * (t155 * t54 + t60 * t272))) * t239) * t233, (((-t12 * t155 + t152 * t24) * t123 + t111 * (t12 * t152 + t155 * t24)) * t153 - ((-t117 * t307 + t155 * t36) * t123 - (t117 * t301 + t152 * t36) * t111) * t239) * t233, (((-t152 * t18 + t155 * t6) * t123 - t111 * (t152 * t6 + t155 * t18)) * t153 + ((-t129 * t307 + t155 * t30) * t123 - (t129 * t301 + t152 * t30) * t111) * t239) * t233; t140 * (t222 * t153 + t85 * t239) * t233, ((t223 * t128 + t222 * t296) * t153 - (t128 * t91 - t85 * t296) * t239) * t233, ((t223 * t116 - t222 * t290) * t153 + (-t116 * t91 - t85 * t290) * t239) * t233, ((((-t53 * t315 + t59 * t79) * t122 + t110 * (-t53 * t316 - t59 * t80)) * t128 + t116 * ((-t59 * t315 - t53 * t79) * t122 - (t59 * t316 - t53 * t80) * t110)) * t153 + (((t152 * t59 + t53 * t274) * t122 + t110 * (t155 * t59 - t53 * t275)) * t128 + t116 * ((-t152 * t53 + t59 * t274) * t122 - t110 * (t155 * t53 + t59 * t275))) * t239) * t233, (((-t11 * t155 + t152 * t23) * t122 + t110 * (t11 * t152 + t155 * t23)) * t153 - ((-t116 * t308 + t155 * t35) * t122 - (t116 * t302 + t152 * t35) * t110) * t239) * t233, (((-t152 * t17 + t155 * t5) * t122 - t110 * (t152 * t5 + t155 * t17)) * t153 + ((-t128 * t308 + t155 * t29) * t122 - (t128 * t302 + t152 * t29) * t110) * t239) * t233; t139 * (t224 * t153 + t84 * t239) * t233, ((t225 * t127 + t224 * t297) * t153 - (t127 * t90 - t84 * t297) * t239) * t233, ((t225 * t115 - t224 * t291) * t153 + (-t115 * t90 - t84 * t291) * t239) * t233, ((((-t52 * t317 + t58 * t79) * t121 + t109 * (-t52 * t318 - t58 * t80)) * t127 + t115 * ((-t58 * t317 - t79 * t52) * t121 - (t58 * t318 - t80 * t52) * t109)) * t153 + (((t152 * t58 + t52 * t277) * t121 + t109 * (t155 * t58 - t52 * t278)) * t127 + t115 * ((-t152 * t52 + t58 * t277) * t121 - t109 * (t52 * t155 + t58 * t278))) * t239) * t233, (((-t10 * t155 + t152 * t22) * t121 + t109 * (t10 * t152 + t155 * t22)) * t153 - ((-t115 * t309 + t155 * t34) * t121 - (t115 * t303 + t152 * t34) * t109) * t239) * t233, (((-t152 * t16 + t155 * t4) * t121 - t109 * (t152 * t4 + t155 * t16)) * t153 + ((-t127 * t309 + t155 * t28) * t121 - (t127 * t303 + t152 * t28) * t109) * t239) * t233; t138 * (t226 * t153 + t83 * t239) * t233, ((t227 * t126 + t226 * t298) * t153 - (t126 * t89 - t83 * t298) * t239) * t233, ((t227 * t114 - t226 * t292) * t153 + (-t114 * t89 - t83 * t292) * t239) * t233, ((((-t51 * t319 + t57 * t79) * t120 + t108 * (-t51 * t320 - t57 * t80)) * t126 + t114 * ((-t57 * t319 - t51 * t79) * t120 - (t57 * t320 - t80 * t51) * t108)) * t153 + (((t152 * t57 + t51 * t280) * t120 + t108 * (t155 * t57 - t51 * t281)) * t126 + t114 * ((-t152 * t51 + t57 * t280) * t120 - t108 * (t155 * t51 + t57 * t281))) * t239) * t233, (((t152 * t21 - t155 * t9) * t120 + t108 * (t152 * t9 + t155 * t21)) * t153 - ((-t114 * t310 + t155 * t33) * t120 - (t114 * t304 + t152 * t33) * t108) * t239) * t233, (((-t15 * t152 + t155 * t3) * t120 - t108 * (t15 * t155 + t152 * t3)) * t153 + ((-t126 * t310 + t155 * t27) * t120 - (t126 * t304 + t152 * t27) * t108) * t239) * t233; t137 * (t228 * t153 + t82 * t239) * t233, ((t229 * t125 + t228 * t299) * t153 - (t125 * t88 - t82 * t299) * t239) * t233, ((t229 * t113 - t228 * t293) * t153 + (-t113 * t88 - t82 * t293) * t239) * t233, ((((-t50 * t321 + t56 * t79) * t119 + t107 * (-t50 * t322 - t56 * t80)) * t125 + t113 * ((-t56 * t321 - t50 * t79) * t119 - (t56 * t322 - t50 * t80) * t107)) * t153 + (((t152 * t56 + t50 * t283) * t119 + t107 * (t155 * t56 - t50 * t284)) * t125 + t113 * ((-t152 * t50 + t56 * t283) * t119 - t107 * (t155 * t50 + t56 * t284))) * t239) * t233, (((t152 * t20 - t155 * t8) * t119 + t107 * (t152 * t8 + t155 * t20)) * t153 - ((-t113 * t311 + t155 * t32) * t119 - (t113 * t305 + t152 * t32) * t107) * t239) * t233, (((-t14 * t152 + t155 * t2) * t119 - t107 * (t14 * t155 + t152 * t2)) * t153 + ((-t125 * t311 + t155 * t26) * t119 - (t125 * t305 + t152 * t26) * t107) * t239) * t233; t136 * (t230 * t153 + t81 * t239) * t233, ((t231 * t124 + t230 * t300) * t153 - (t124 * t87 - t81 * t300) * t239) * t233, ((t231 * t112 - t230 * t294) * t153 + (-t112 * t87 - t81 * t294) * t239) * t233, ((((-t49 * t323 + t55 * t79) * t118 + t106 * (-t49 * t324 - t55 * t80)) * t124 + t112 * ((-t55 * t323 - t49 * t79) * t118 - (t55 * t324 - t49 * t80) * t106)) * t153 + (((t152 * t55 + t49 * t286) * t118 + t106 * (t155 * t55 - t49 * t287)) * t124 + t112 * ((-t152 * t49 + t55 * t286) * t118 - t106 * (t155 * t49 + t55 * t287))) * t239) * t233, (((t152 * t19 - t155 * t7) * t118 + t106 * (t152 * t7 + t155 * t19)) * t153 - ((-t112 * t312 + t155 * t31) * t118 - (t112 * t306 + t152 * t31) * t106) * t239) * t233, (((t1 * t155 - t13 * t152) * t118 - t106 * (t1 * t152 + t13 * t155)) * t153 + ((-t124 * t312 + t155 * t25) * t118 - (t124 * t306 + t152 * t25) * t106) * t239) * t233;];
Jinv  = t61;
