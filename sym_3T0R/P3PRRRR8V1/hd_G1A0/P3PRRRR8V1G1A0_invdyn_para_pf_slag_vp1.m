% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRRR8V1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 16:50
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRRR8V1G1A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 16:49:37
% EndTime: 2020-08-06 16:49:43
% DurationCPUTime: 5.39s
% Computational Cost: add. (15780->381), mult. (41692->747), div. (4053->9), fcn. (48402->28), ass. (0->315)
t186 = cos(qJ(3,2));
t165 = 0.1e1 / t186;
t187 = cos(qJ(2,2));
t181 = sin(qJ(2,2));
t265 = t181 * t186;
t135 = pkin(2) * t265 - t187 * pkin(5);
t169 = sin(pkin(3));
t171 = cos(pkin(3));
t180 = sin(qJ(3,2));
t229 = 0.1e1 / (t171 * t180 * pkin(2) + t135 * t169);
t293 = t229 * t165;
t324 = g(3) * t171;
t182 = sin(qJ(3,1));
t188 = cos(qJ(3,1));
t220 = rSges(3,1) * t188 - rSges(3,2) * t182;
t354 = t220 * m(3);
t221 = rSges(3,1) * t186 - rSges(3,2) * t180;
t353 = t221 * m(3);
t178 = sin(qJ(3,3));
t184 = cos(qJ(3,3));
t222 = rSges(3,1) * t184 - rSges(3,2) * t178;
t352 = t222 * m(3);
t174 = legFrame(1,3);
t157 = sin(t174);
t160 = cos(t174);
t130 = -t157 * g(1) + t160 * g(2);
t133 = t160 * g(1) + t157 * g(2);
t168 = sin(pkin(6));
t170 = cos(pkin(6));
t212 = t130 * t170 - t133 * t168;
t325 = g(3) * t169;
t351 = t212 * t171 + t325;
t173 = legFrame(2,3);
t156 = sin(t173);
t159 = cos(t173);
t129 = -t156 * g(1) + t159 * g(2);
t132 = t159 * g(1) + t156 * g(2);
t214 = t129 * t170 - t132 * t168;
t350 = t214 * t171 + t325;
t172 = legFrame(3,3);
t155 = sin(t172);
t158 = cos(t172);
t128 = -t155 * g(1) + t158 * g(2);
t131 = t158 * g(1) + t155 * g(2);
t216 = t128 * t170 - t131 * t168;
t349 = t216 * t171 + t325;
t179 = sin(qJ(2,3));
t268 = t179 * t184;
t238 = t169 * t268;
t276 = t171 * t178;
t185 = cos(qJ(2,3));
t281 = t169 * t185;
t348 = 0.1e1 / (-pkin(5) * t281 + (t238 + t276) * pkin(2));
t183 = sin(qJ(2,1));
t262 = t183 * t188;
t236 = t169 * t262;
t273 = t171 * t182;
t189 = cos(qJ(2,1));
t277 = t169 * t189;
t347 = 0.1e1 / (-pkin(5) * t277 + (t236 + t273) * pkin(2));
t154 = m(3) * rSges(3,1) * rSges(3,2) - Icges(3,4);
t346 = 0.2e1 * t154;
t193 = m(2) * rSges(2,1);
t345 = m(3) * rSges(3,3);
t140 = t178 * rSges(3,1) + t184 * rSges(3,2);
t94 = -t169 * t179 * t140 + t222 * t171;
t344 = m(3) * t94;
t141 = t180 * rSges(3,1) + t186 * rSges(3,2);
t95 = -t169 * t181 * t141 + t221 * t171;
t343 = m(3) * t95;
t142 = t182 * rSges(3,1) + t188 * rSges(3,2);
t96 = -t169 * t183 * t142 + t220 * t171;
t342 = m(3) * t96;
t194 = xDP(2);
t195 = xDP(1);
t203 = 0.1e1 / pkin(2);
t163 = 0.1e1 / t184;
t294 = t348 * t163;
t241 = t203 * t294;
t116 = t170 * t155 + t158 * t168;
t260 = t185 * t116;
t113 = -t168 * t155 + t158 * t170;
t261 = t185 * t113;
t269 = t179 * t116;
t270 = t179 * t113;
t329 = pkin(2) * t184;
t70 = -(t171 * t261 - t269) * t329 - pkin(5) * (t171 * t270 + t260);
t71 = -(t171 * t260 + t270) * t329 - (t171 * t269 - t261) * pkin(5);
t55 = (t194 * t71 + t195 * t70) * t241;
t341 = pkin(2) * t55;
t240 = t203 * t293;
t117 = t170 * t156 + t159 * t168;
t258 = t187 * t117;
t114 = -t168 * t156 + t159 * t170;
t259 = t187 * t114;
t266 = t181 * t117;
t267 = t181 * t114;
t328 = pkin(2) * t186;
t72 = -(t171 * t259 - t266) * t328 - pkin(5) * (t171 * t267 + t258);
t73 = -(t171 * t258 + t267) * t328 - (t171 * t266 - t259) * pkin(5);
t56 = (t194 * t73 + t195 * t72) * t240;
t340 = pkin(2) * t56;
t167 = 0.1e1 / t188;
t292 = t347 * t167;
t239 = t203 * t292;
t118 = t170 * t157 + t160 * t168;
t256 = t189 * t118;
t115 = -t168 * t157 + t160 * t170;
t257 = t189 * t115;
t263 = t183 * t118;
t264 = t183 * t115;
t327 = pkin(2) * t188;
t74 = -(t171 * t257 - t263) * t327 - pkin(5) * (t171 * t264 + t256);
t75 = -(t171 * t256 + t264) * t327 - (t171 * t263 - t257) * pkin(5);
t57 = (t194 * t75 + t195 * t74) * t239;
t339 = pkin(2) * t57;
t134 = pkin(2) * t268 - t185 * pkin(5);
t107 = 0.1e1 / (pkin(2) * t276 + t134 * t169);
t275 = t171 * t179;
t119 = t168 * t185 + t170 * t275;
t122 = -t168 * t275 + t170 * t185;
t282 = t169 * t184;
t76 = (-t119 * t158 - t155 * t122) * t178 - t113 * t282;
t79 = (-t155 * t119 + t122 * t158) * t178 - t116 * t282;
t61 = (t194 * t79 + t195 * t76) * t163 * t107;
t338 = pkin(5) * t61;
t274 = t171 * t181;
t120 = t168 * t187 + t170 * t274;
t123 = -t168 * t274 + t170 * t187;
t280 = t169 * t186;
t77 = (-t120 * t159 - t156 * t123) * t180 - t114 * t280;
t80 = (-t156 * t120 + t123 * t159) * t180 - t117 * t280;
t62 = (t194 * t80 + t195 * t77) * t293;
t337 = pkin(5) * t62;
t136 = pkin(2) * t262 - t189 * pkin(5);
t109 = 0.1e1 / (pkin(2) * t273 + t136 * t169);
t272 = t171 * t183;
t121 = t168 * t189 + t170 * t272;
t124 = -t168 * t272 + t170 * t189;
t278 = t169 * t188;
t78 = (-t121 * t160 - t157 * t124) * t182 - t115 * t278;
t81 = (-t157 * t121 + t124 * t160) * t182 - t118 * t278;
t63 = (t194 * t81 + t195 * t78) * t167 * t109;
t336 = pkin(5) * t63;
t199 = rSges(3,2) ^ 2;
t200 = rSges(3,1) ^ 2;
t144 = (-t199 + t200) * m(3) - Icges(3,1) + Icges(3,2);
t335 = t144 / 0.2e1;
t152 = rSges(3,2) * t345 - Icges(3,6);
t334 = -t152 / 0.4e1;
t153 = rSges(3,1) * t345 - Icges(3,5);
t333 = t153 / 0.4e1;
t162 = t184 ^ 2;
t332 = pkin(2) * t162;
t164 = t186 ^ 2;
t331 = pkin(2) * t164;
t166 = t188 ^ 2;
t330 = pkin(2) * t166;
t161 = m(1) + m(2) + m(3);
t326 = g(3) * t161;
t201 = pkin(5) ^ 2;
t202 = pkin(2) ^ 2;
t296 = t180 * t56;
t251 = pkin(2) * t296;
t323 = (-pkin(5) * t251 + (t164 * t202 + t201) * t62) * t62;
t249 = t178 * t338;
t297 = t178 * t55;
t252 = pkin(2) * t297;
t25 = -pkin(5) * t252 + (t162 * t202 + t201) * t61;
t271 = t171 * t203;
t285 = t169 * t178;
t314 = t348 * t61;
t16 = (t25 * t271 * t314 + (-t55 * t134 * t285 + t171 * (t55 * t332 - t249)) * t107 * t55) * t163;
t322 = t94 * t16;
t248 = t180 * t337;
t284 = t169 * t180;
t17 = (t271 * t323 + (-t56 * t135 * t284 + t171 * (t56 * t331 - t248)) * t56) * t293;
t321 = t95 * t17;
t247 = t182 * t336;
t295 = t182 * t57;
t250 = pkin(2) * t295;
t27 = -pkin(5) * t250 + (t166 * t202 + t201) * t63;
t283 = t169 * t182;
t313 = t347 * t63;
t18 = (t27 * t271 * t313 + (-t57 * t136 * t283 + t171 * (t57 * t330 - t247)) * t109 * t57) * t167;
t320 = t96 * t18;
t151 = m(2) * rSges(2,2) - t345;
t28 = t249 - t341;
t22 = (t25 * t61 - t28 * t341) * t348;
t253 = 0.2e1 * m(3);
t52 = t55 ^ 2;
t58 = t61 ^ 2;
t319 = ((-t58 * t193 - (t58 + t52) * t352) * t179 - (t140 * t55 * t253 + t61 * t151) * t185 * t61) * t169 + t161 * t22;
t29 = t248 - t340;
t23 = (t29 * t340 - t323) * t229;
t53 = t56 ^ 2;
t59 = t62 ^ 2;
t318 = ((-t59 * t193 - (t59 + t53) * t353) * t181 - (t141 * t56 * t253 + t62 * t151) * t187 * t62) * t169 - t161 * t23;
t30 = -t247 + t339;
t24 = (t27 * t63 + t30 * t339) * t347;
t54 = t57 ^ 2;
t60 = t63 ^ 2;
t317 = ((-t60 * t193 - (t60 + t54) * t354) * t183 - (t142 * t57 * t253 + t63 * t151) * t189 * t63) * t169 + t161 * t24;
t315 = rSges(3,2) * t169;
t312 = t140 * t52;
t311 = t141 * t53;
t310 = t142 * t54;
t309 = t163 * t76;
t308 = t163 * t79;
t196 = 0.2e1 * qJ(3,3);
t254 = t199 + t200;
t207 = Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1 + (0.2e1 * rSges(3,3) ^ 2 + t254) * m(3) / 0.2e1;
t82 = cos(t196) * t335 - t154 * sin(t196) + t207;
t307 = t163 * t82;
t306 = t165 * t77;
t305 = t165 * t80;
t197 = 0.2e1 * qJ(3,2);
t83 = cos(t197) * t335 - t154 * sin(t197) + t207;
t304 = t165 * t83;
t303 = t167 * t78;
t302 = t167 * t81;
t198 = 0.2e1 * qJ(3,1);
t84 = cos(t198) * t335 - t154 * sin(t198) + t207;
t301 = t167 * t84;
t125 = t193 + t352;
t97 = t125 * t185 - t179 * t151;
t300 = t169 * t97;
t126 = t193 + t353;
t98 = t126 * t187 - t181 * t151;
t299 = t169 * t98;
t127 = t193 + t354;
t99 = t127 * t189 - t183 * t151;
t298 = t169 * t99;
t110 = -t152 * t184 - t153 * t178;
t291 = t110 * t163;
t111 = -t152 * t186 - t153 * t180;
t290 = t111 * t165;
t112 = -t152 * t188 - t153 * t182;
t289 = t112 * t167;
t288 = t144 * t178;
t287 = t144 * t180;
t286 = t144 * t182;
t279 = t169 * t187;
t255 = rSges(3,2) * t324;
t228 = t241 * t344;
t244 = t163 * t300;
t137 = pkin(5) * t179 + t185 * t329;
t210 = pkin(2) * t285 - t134 * t171;
t85 = t137 * t170 + t210 * t168;
t88 = t168 * t137 - t210 * t170;
t64 = -t155 * t88 + t85 * t158;
t40 = t70 * t228 + (t161 * t64 + t76 * t244) * t107;
t227 = t240 * t343;
t243 = t165 * t299;
t138 = pkin(5) * t181 + t187 * t328;
t209 = pkin(2) * t284 - t135 * t171;
t86 = t138 * t170 + t209 * t168;
t89 = t168 * t138 - t209 * t170;
t65 = -t156 * t89 + t86 * t159;
t41 = t72 * t227 + (t161 * t65 + t77 * t243) * t229;
t226 = t239 * t342;
t242 = t167 * t298;
t139 = pkin(5) * t183 + t189 * t327;
t208 = pkin(2) * t283 - t136 * t171;
t87 = t139 * t170 + t208 * t168;
t90 = t168 * t139 - t208 * t170;
t66 = -t157 * t90 + t87 * t160;
t42 = t74 * t226 + (t161 * t66 + t78 * t242) * t109;
t246 = t42 + t41 + t40;
t67 = t85 * t155 + t88 * t158;
t43 = t71 * t228 + (t161 * t67 + t79 * t244) * t107;
t68 = t86 * t156 + t89 * t159;
t44 = t73 * t227 + (t161 * t68 + t80 * t243) * t229;
t69 = t87 * t157 + t90 * t160;
t45 = t75 * t226 + (t161 * t69 + t81 * t242) * t109;
t245 = t45 + t44 + t43;
t237 = t169 * t265;
t235 = t110 * t241;
t149 = t254 * m(3) + Icges(3,3);
t234 = t149 * t241;
t233 = t111 * t240;
t232 = t149 * t240;
t231 = t112 * t239;
t230 = t149 * t239;
t176 = xDDP(2);
t177 = xDDP(1);
t219 = t176 * t71 + t177 * t70;
t218 = t176 * t73 + t177 * t72;
t217 = t176 * t75 + t177 * t74;
t215 = t168 * t128 + t131 * t170;
t213 = t168 * t129 + t132 * t170;
t211 = t168 * t130 + t133 * t170;
t206 = t349 * t179 + t215 * t185;
t205 = t350 * t181 + t213 * t187;
t204 = t351 * t183 + t211 * t189;
t175 = xDDP(3);
t39 = t75 * t231 + (t69 * t298 + t81 * t301) * t109;
t38 = t73 * t233 + (t68 * t299 + t80 * t304) * t229;
t37 = t71 * t235 + (t67 * t300 + t79 * t307) * t107;
t36 = t74 * t231 + (t66 * t298 + t78 * t301) * t109;
t35 = t72 * t233 + (t65 * t299 + t77 * t304) * t229;
t34 = t70 * t235 + (t64 * t300 + t76 * t307) * t107;
t12 = (((t171 * t57 + t63 * t277) * t330 - (t250 - t336) * t236 - t171 * t30) * t313 + (t57 * t277 + (t166 * t171 - t182 * t236 - t171) * t63) * t347 * t339) * t167;
t11 = (((t171 * t56 + t62 * t279) * t331 - (t251 - t337) * t237 + t171 * t29) * t62 + (t56 * t279 + (t164 * t171 - t180 * t237 - t171) * t62) * t340) * t293;
t10 = (((t171 * t55 + t61 * t281) * t332 - (t252 - t338) * t238 + t171 * t28) * t314 + (t55 * t281 + (t162 * t171 - t178 * t238 - t171) * t61) * t348 * t341) * t163;
t9 = t24 * t342 - t112 * t12 - t149 * t18 + (t166 * t346 + t188 * t286 - t154) * t60 + m(3) * (((t212 * t169 - t324) * rSges(3,1) + t204 * rSges(3,2)) * t188 + (t204 * rSges(3,1) - t212 * t315 + t255) * t182);
t8 = -t23 * t343 - t111 * t11 - t149 * t17 + (t164 * t346 + t186 * t287 - t154) * t59 + m(3) * (((t214 * t169 - t324) * rSges(3,1) + t205 * rSges(3,2)) * t186 + (t205 * rSges(3,1) - t214 * t315 + t255) * t180);
t7 = t22 * t344 - t110 * t10 - t149 * t16 + (t162 * t346 + t184 * t288 - t154) * t58 + m(3) * (((t216 * t169 - t324) * rSges(3,1) + t206 * rSges(3,2)) * t184 + (t206 * rSges(3,1) - t216 * t315 + t255) * t178);
t6 = t24 * t298 - t84 * t12 - t112 * t18 - 0.4e1 * t57 * ((t63 * t286 / 0.2e1 + t57 * t333) * t188 + t295 * t334 + (t166 - 0.1e1 / 0.2e1) * t63 * t154) + (-t127 * t351 + t211 * t151) * t189 - t183 * (-t211 * t127 - t151 * t351);
t5 = -t23 * t299 - t83 * t11 - t111 * t17 - 0.4e1 * t56 * ((t62 * t287 / 0.2e1 + t56 * t333) * t186 + t296 * t334 + (t164 - 0.1e1 / 0.2e1) * t62 * t154) + (-t126 * t350 + t213 * t151) * t187 - t181 * (-t213 * t126 - t151 * t350);
t4 = t22 * t300 - t82 * t10 - t110 * t16 - 0.4e1 * t55 * ((t61 * t288 / 0.2e1 + t55 * t333) * t184 + t297 * t334 + (t162 - 0.1e1 / 0.2e1) * t61 * t154) + (-t125 * t349 + t215 * t151) * t185 - t179 * (-t215 * t125 - t151 * t349);
t3 = -t12 * t298 - t326 + (-t171 * t310 - t320) * m(3) + t317;
t2 = -t11 * t299 - t326 + (-t171 * t311 - t321) * m(3) + t318;
t1 = -t10 * t300 - t326 + (-t171 * t312 - t322) * m(3) + t319;
t13 = [(-g(1) + t177) * m(4) + t246 * t175 + ((t36 * t303 + t42 * t66) * t177 + (t36 * t302 + t42 * t69) * t176 + t66 * t3 + t6 * t303) * t109 + ((t35 * t306 + t41 * t65) * t177 + (t35 * t305 + t41 * t68) * t176 + t65 * t2 + t5 * t306) * t229 + ((t74 * t9 + t217 * (t74 * t230 + (t78 * t289 + t66 * t342) * t109)) * t292 + (t72 * t8 + t218 * (t72 * t232 + (t77 * t290 + t65 * t343) * t229)) * t293 + (t7 * t70 + t219 * (t70 * t234 + (t76 * t291 + t64 * t344) * t107)) * t294) * t203 + ((t34 * t309 + t40 * t64) * t177 + (t34 * t308 + t40 * t67) * t176 + t64 * t1 + t4 * t309) * t107; (-g(2) + t176) * m(4) + t245 * t175 + ((t39 * t303 + t45 * t66) * t177 + (t39 * t302 + t45 * t69) * t176 + t69 * t3 + t6 * t302) * t109 + ((t38 * t306 + t44 * t65) * t177 + (t38 * t305 + t44 * t68) * t176 + t68 * t2 + t5 * t305) * t229 + ((t75 * t9 + t217 * (t75 * t230 + (t81 * t289 + t69 * t342) * t109)) * t292 + (t73 * t8 + t218 * (t73 * t232 + (t80 * t290 + t68 * t343) * t229)) * t293 + (t7 * t71 + t219 * (t71 * t234 + (t79 * t291 + t67 * t344) * t107)) * t294) * t203 + ((t37 * t309 + t43 * t64) * t177 + (t37 * t308 + t43 * t67) * t176 + t67 * t1 + t4 * t308) * t107; t246 * t177 + t245 * t176 + (-t97 * t10 - t98 * t11 - t99 * t12) * t169 + (-t322 - t321 - t320 + (-t310 - t311 - t312) * t171) * m(3) + t317 + t318 + t319 + (-t175 + g(3)) * (-0.3e1 * t161 - m(4));];
tauX  = t13;
