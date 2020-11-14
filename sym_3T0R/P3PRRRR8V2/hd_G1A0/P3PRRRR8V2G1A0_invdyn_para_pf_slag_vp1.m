% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRRR8V2G1A0
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2020-08-06 17:36
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRRR8V2G1A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:35:08
% EndTime: 2020-08-06 17:35:15
% DurationCPUTime: 6.26s
% Computational Cost: add. (26460->447), mult. (51760->799), div. (2757->10), fcn. (53592->22), ass. (0->326)
t210 = cos(qJ(2,3));
t221 = pkin(7) + pkin(6);
t171 = t210 * t221;
t204 = sin(qJ(2,3));
t140 = pkin(2) * t204 - t171;
t194 = sin(pkin(4));
t209 = cos(qJ(3,3));
t196 = cos(pkin(4));
t203 = sin(qJ(3,3));
t290 = t196 * t203;
t301 = t194 * t204;
t190 = t209 ^ 2;
t351 = pkin(3) * t190;
t85 = 0.1e1 / ((pkin(3) * t290 + t140 * t194) * t209 + pkin(2) * t290 + t301 * t351);
t212 = cos(qJ(2,2));
t172 = t212 * t221;
t206 = sin(qJ(2,2));
t141 = pkin(2) * t206 - t172;
t211 = cos(qJ(3,2));
t205 = sin(qJ(3,2));
t288 = t196 * t205;
t299 = t194 * t206;
t191 = t211 ^ 2;
t350 = pkin(3) * t191;
t86 = 0.1e1 / ((pkin(3) * t288 + t141 * t194) * t211 + pkin(2) * t288 + t299 * t350);
t214 = cos(qJ(2,1));
t173 = t214 * t221;
t208 = sin(qJ(2,1));
t142 = pkin(2) * t208 - t173;
t213 = cos(qJ(3,1));
t207 = sin(qJ(3,1));
t286 = t196 * t207;
t297 = t194 * t208;
t192 = t213 ^ 2;
t349 = pkin(3) * t192;
t87 = 0.1e1 / ((pkin(3) * t286 + t142 * t194) * t213 + pkin(2) * t286 + t297 * t349);
t240 = rSges(3,1) * t213 - rSges(3,2) * t207;
t379 = m(3) * t240;
t241 = rSges(3,1) * t211 - rSges(3,2) * t205;
t378 = m(3) * t241;
t242 = rSges(3,1) * t209 - rSges(3,2) * t203;
t377 = m(3) * t242;
t346 = g(3) * t196;
t278 = t208 * t221;
t145 = pkin(2) * t214 + t278;
t193 = sin(pkin(8));
t195 = cos(pkin(8));
t298 = t194 * t207;
t231 = pkin(3) * t298 - t142 * t196;
t373 = t145 * t195 + t193 * t231;
t280 = t206 * t221;
t144 = pkin(2) * t212 + t280;
t300 = t194 * t205;
t232 = pkin(3) * t300 - t141 * t196;
t372 = t144 * t195 + t193 * t232;
t282 = t204 * t221;
t143 = pkin(2) * t210 + t282;
t302 = t194 * t203;
t233 = pkin(3) * t302 - t140 * t196;
t371 = t143 * t195 + t193 * t233;
t199 = legFrame(1,3);
t180 = sin(t199);
t183 = cos(t199);
t132 = -t180 * g(1) + t183 * g(2);
t135 = t183 * g(1) + t180 * g(2);
t235 = t132 * t195 - t135 * t193;
t347 = g(3) * t194;
t370 = t235 * t196 + t347;
t198 = legFrame(2,3);
t179 = sin(t198);
t182 = cos(t198);
t131 = -t179 * g(1) + t182 * g(2);
t134 = t182 * g(1) + t179 * g(2);
t237 = t131 * t195 - t134 * t193;
t369 = t237 * t196 + t347;
t197 = legFrame(3,3);
t178 = sin(t197);
t181 = cos(t197);
t130 = -t178 * g(1) + t181 * g(2);
t133 = t181 * g(1) + t178 * g(2);
t239 = t130 * t195 - t133 * t193;
t368 = t239 * t196 + t347;
t362 = m(3) * rSges(3,1);
t269 = rSges(3,2) * t362;
t174 = -Icges(3,4) + t269;
t222 = pkin(2) * m(3);
t265 = t222 / 0.2e1;
t250 = rSges(3,1) * t265;
t367 = t174 * t190 + t203 * t250;
t366 = t174 * t191 + t205 * t250;
t365 = t174 * t192 + t207 * t250;
t364 = 0.2e1 * pkin(2);
t363 = -0.2e1 * t174;
t219 = xDP(2);
t220 = xDP(1);
t226 = 0.1e1 / pkin(3);
t168 = t209 * pkin(3) + pkin(2);
t127 = t204 * t168 - t171;
t136 = t168 * t290;
t296 = t194 * t209;
t91 = 0.1e1 / (t127 * t296 + t136);
t318 = t226 * t91;
t109 = -t193 * t178 + t181 * t195;
t112 = t195 * t178 + t181 * t193;
t312 = (t168 * t210 + t282) * t196;
t76 = -t127 * t109 - t112 * t312;
t79 = -t109 * t312 + t127 * t112;
t61 = (t219 * t76 + t220 * t79) * t318;
t361 = pkin(3) * t61;
t169 = t211 * pkin(3) + pkin(2);
t128 = t206 * t169 - t172;
t137 = t169 * t288;
t294 = t194 * t211;
t92 = 0.1e1 / (t128 * t294 + t137);
t317 = t226 * t92;
t110 = -t193 * t179 + t182 * t195;
t113 = t195 * t179 + t182 * t193;
t311 = (t169 * t212 + t280) * t196;
t77 = -t128 * t110 - t113 * t311;
t80 = -t110 * t311 + t128 * t113;
t62 = (t219 * t77 + t220 * t80) * t317;
t360 = pkin(3) * t62;
t170 = t213 * pkin(3) + pkin(2);
t129 = t208 * t170 - t173;
t138 = t170 * t286;
t292 = t194 * t213;
t93 = 0.1e1 / (t129 * t292 + t138);
t316 = t226 * t93;
t111 = -t193 * t180 + t183 * t195;
t114 = t195 * t180 + t183 * t193;
t310 = (t170 * t214 + t278) * t196;
t78 = -t129 * t111 - t114 * t310;
t81 = -t111 * t310 + t129 * t114;
t63 = (t219 * t78 + t220 * t81) * t316;
t359 = pkin(3) * t63;
t223 = rSges(3,2) ^ 2;
t224 = rSges(3,1) ^ 2;
t155 = (-t223 + t224) * m(3) + Icges(3,2) - Icges(3,1);
t358 = t155 / 0.2e1;
t215 = pkin(6) + rSges(3,3);
t345 = t215 * m(3);
t159 = rSges(3,2) * t345 - Icges(3,6);
t357 = -t159 / 0.4e1;
t160 = rSges(3,1) * t345 - Icges(3,5);
t356 = t160 / 0.4e1;
t355 = -t174 / 0.2e1;
t152 = t203 * rSges(3,1) + t209 * rSges(3,2);
t100 = -t152 * t301 + t196 * t242;
t354 = m(3) * t100;
t153 = t205 * rSges(3,1) + t211 * rSges(3,2);
t101 = -t153 * t299 + t196 * t241;
t353 = m(3) * t101;
t154 = t207 * rSges(3,1) + t213 * rSges(3,2);
t102 = -t154 * t297 + t196 * t240;
t352 = m(3) * t102;
t189 = m(1) + m(2) + m(3);
t348 = g(3) * t189;
t227 = pkin(2) ^ 2;
t175 = t221 ^ 2 + t227;
t225 = pkin(3) ^ 2;
t268 = t203 * t361;
t276 = pkin(3) * t364;
t289 = t196 * t204;
t70 = -t109 * t296 - (t109 * t289 + t210 * t112) * t203;
t73 = -t112 * t296 - (-t210 * t109 + t112 * t289) * t203;
t55 = (t219 * t73 + t220 * t70) * t85;
t344 = (-t221 * t268 + (t190 * t225 + t209 * t276 + t175) * t55) * t55;
t267 = t205 * t360;
t287 = t196 * t206;
t71 = -t110 * t294 - (t110 * t287 + t212 * t113) * t205;
t74 = -t113 * t294 - (-t212 * t110 + t113 * t287) * t205;
t56 = (t219 * t74 + t220 * t71) * t86;
t343 = (-t221 * t267 + (t191 * t225 + t211 * t276 + t175) * t56) * t56;
t266 = t207 * t359;
t285 = t196 * t208;
t72 = -t111 * t292 - (t111 * t285 + t214 * t114) * t207;
t75 = -t114 * t292 - (-t214 * t111 + t114 * t285) * t207;
t57 = (t219 * t75 + t220 * t72) * t87;
t342 = (-t221 * t266 + (t192 * t225 + t213 * t276 + t175) * t57) * t57;
t341 = t76 * t91;
t340 = t77 * t92;
t339 = t78 * t93;
t338 = t79 * t91;
t337 = t80 * t92;
t336 = t81 * t93;
t157 = m(2) * rSges(2,2) - t345;
t167 = m(2) * rSges(2,1) + t222;
t315 = t55 * t221;
t253 = t203 * t315;
t46 = t253 - t361;
t19 = (-t209 * t344 - (pkin(2) * t61 - t46 * t209) * t361) * t85;
t270 = 0.2e1 * m(3);
t321 = t210 * t55;
t52 = t55 ^ 2;
t58 = t61 ^ 2;
t334 = -t189 * t19 + ((-t52 * t167 - (t52 + t58) * t377) * t204 - (t152 * t61 * t270 + t157 * t55) * t321) * t194;
t314 = t56 * t221;
t252 = t205 * t314;
t47 = t252 - t360;
t20 = (-t211 * t343 - (pkin(2) * t62 - t47 * t211) * t360) * t86;
t320 = t212 * t56;
t53 = t56 ^ 2;
t59 = t62 ^ 2;
t333 = -t189 * t20 + ((-t53 * t167 - (t53 + t59) * t378) * t206 - (t153 * t62 * t270 + t157 * t56) * t320) * t194;
t313 = t57 * t221;
t251 = t207 * t313;
t48 = t251 - t359;
t21 = (-t213 * t342 - (pkin(2) * t63 - t48 * t213) * t359) * t87;
t319 = t214 * t57;
t54 = t57 ^ 2;
t60 = t63 ^ 2;
t332 = -t189 * t21 + ((-t54 * t167 - (t54 + t60) * t379) * t208 - (t154 * t63 * t270 + t157 * t57) * t319) * t194;
t331 = rSges(3,2) * t194;
t283 = t204 * t209;
t106 = pkin(3) * t283 + t140;
t284 = t196 * t226;
t13 = t85 * t284 * t344 + (-t196 * t253 + (-t106 * t302 + (pkin(2) * t209 + t351) * t196) * t61) / (t106 * t296 + t136) * t61;
t330 = t100 * t13;
t281 = t206 * t211;
t107 = pkin(3) * t281 + t141;
t14 = t86 * t284 * t343 + (-t196 * t252 + (-t107 * t300 + (pkin(2) * t211 + t350) * t196) * t62) / (t107 * t294 + t137) * t62;
t329 = t101 * t14;
t279 = t208 * t213;
t108 = pkin(3) * t279 + t142;
t15 = t87 * t284 * t342 + (-t196 * t251 + (-t108 * t298 + (pkin(2) * t213 + t349) * t196) * t63) / (t108 * t292 + t138) * t63;
t328 = t102 * t15;
t327 = t152 * t58;
t326 = t153 * t59;
t325 = t154 * t60;
t115 = t167 + t377;
t97 = t115 * t210 - t157 * t204;
t324 = t194 * t97;
t116 = t167 + t378;
t98 = t116 * t212 - t157 * t206;
t323 = t194 * t98;
t117 = t167 + t379;
t99 = t117 * t214 - t157 * t208;
t322 = t194 * t99;
t156 = (t223 + t224) * m(3) + Icges(3,3);
t306 = t156 * t226;
t295 = t194 * t210;
t293 = t194 * t212;
t291 = t194 * t214;
t277 = rSges(3,2) * t346;
t272 = -t269 / 0.2e1 + Icges(3,4) / 0.2e1;
t271 = -0.2e1 * rSges(3,2) * pkin(2);
t249 = t318 * t354;
t118 = t193 * t289 - t195 * t210;
t121 = t193 * t210 + t195 * t289;
t262 = pkin(2) * t302;
t88 = t193 * t143 - t195 * t233;
t64 = -(t118 * t181 + t178 * t121) * t351 + (-t88 * t178 + t371 * t181) * t209 + t112 * t262;
t34 = t79 * t249 + (t189 * t64 + t324 * t70) * t85;
t248 = t317 * t353;
t119 = t193 * t287 - t195 * t212;
t122 = t193 * t212 + t195 * t287;
t261 = pkin(2) * t300;
t89 = t193 * t144 - t195 * t232;
t65 = -(t119 * t182 + t179 * t122) * t350 + (-t89 * t179 + t372 * t182) * t211 + t113 * t261;
t35 = t80 * t248 + (t189 * t65 + t323 * t71) * t86;
t247 = t316 * t352;
t120 = t193 * t285 - t195 * t214;
t123 = t193 * t214 + t195 * t285;
t260 = pkin(2) * t298;
t90 = t193 * t145 - t195 * t231;
t66 = -(t120 * t183 + t180 * t123) * t349 + (-t90 * t180 + t373 * t183) * t213 + t114 * t260;
t36 = t81 * t247 + (t189 * t66 + t322 * t72) * t87;
t264 = t36 + t35 + t34;
t67 = (-t178 * t118 + t121 * t181) * t351 + (t371 * t178 + t88 * t181) * t209 - t109 * t262;
t37 = t76 * t249 + (t189 * t67 + t324 * t73) * t85;
t68 = (-t179 * t119 + t122 * t182) * t350 + (t372 * t179 + t89 * t182) * t211 - t110 * t261;
t38 = t77 * t248 + (t189 * t68 + t323 * t74) * t86;
t69 = (-t180 * t120 + t123 * t183) * t349 + (t373 * t180 + t90 * t183) * t213 - t111 * t260;
t39 = t78 * t247 + (t189 * t69 + t322 * t75) * t87;
t263 = t39 + t38 + t37;
t103 = -t159 * t209 - t160 * t203;
t259 = t103 * t318;
t104 = -t159 * t211 - t160 * t205;
t258 = t104 * t317;
t105 = -t159 * t213 - t160 * t207;
t257 = t105 * t316;
t256 = t91 * t306;
t255 = t92 * t306;
t254 = t93 * t306;
t177 = rSges(3,2) * t265;
t246 = (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,1) + Icges(2,3);
t238 = t193 * t130 + t133 * t195;
t236 = t193 * t131 + t134 * t195;
t234 = t193 * t132 + t135 * t195;
t230 = t368 * t204 + t238 * t210;
t229 = t369 * t206 + t236 * t212;
t228 = t370 * t208 + t234 * t214;
t202 = xDDP(1);
t201 = xDDP(2);
t200 = xDDP(3);
t185 = t362 * t364;
t158 = t215 ^ 2 + t223 + t227;
t139 = (t224 / 0.2e1 - t223 / 0.2e1) * m(3) - Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1;
t84 = t155 * t192 + (t207 * t363 + t185) * t213 + (t207 * t271 + t158) * m(3) + t246;
t83 = t155 * t191 + (t205 * t363 + t185) * t211 + (t205 * t271 + t158) * m(3) + t246;
t82 = t155 * t190 + (t203 * t363 + t185) * t209 + (t203 * t271 + t158) * m(3) + t246;
t45 = t78 * t254 + (t105 * t75 + t352 * t69) * t87;
t44 = t77 * t255 + (t104 * t74 + t353 * t68) * t86;
t43 = t76 * t256 + (t103 * t73 + t354 * t67) * t85;
t42 = t81 * t254 + (t105 * t72 + t352 * t66) * t87;
t41 = t80 * t255 + (t104 * t71 + t353 * t65) * t86;
t40 = t79 * t256 + (t103 * t70 + t354 * t64) * t85;
t33 = t78 * t257 + (t322 * t69 + t75 * t84) * t87;
t32 = t77 * t258 + (t323 * t68 + t74 * t83) * t86;
t31 = t76 * t259 + (t324 * t67 + t73 * t82) * t85;
t30 = t81 * t257 + (t322 * t66 + t72 * t84) * t87;
t29 = t80 * t258 + (t323 * t65 + t71 * t83) * t86;
t28 = t79 * t259 + (t324 * t64 + t70 * t82) * t85;
t12 = (((t196 * t63 + t291 * t57) * t349 + ((-t266 + t313) * t208 + pkin(2) * t319) * t292 + t196 * t48) * t57 + (t63 * t291 + (t192 * t196 - t279 * t298 - t196) * t57) * t359) * t87;
t11 = (((t196 * t62 + t293 * t56) * t350 + ((-t267 + t314) * t206 + pkin(2) * t320) * t294 + t196 * t47) * t56 + (t62 * t293 + (t191 * t196 - t281 * t300 - t196) * t56) * t360) * t86;
t10 = (((t196 * t61 + t295 * t55) * t351 + ((-t268 + t315) * t204 + pkin(2) * t321) * t296 + t196 * t46) * t55 + (t61 * t295 + (t190 * t196 - t283 * t302 - t196) * t55) * t361) * t85;
t9 = -t21 * t352 - t105 * t12 - t156 * t15 + 0.2e1 * ((t139 * t207 + t177) * t213 + t272 + t365) * t54 + (((t194 * t235 - t346) * rSges(3,1) + t228 * rSges(3,2)) * t213 + t207 * (rSges(3,1) * t228 - t235 * t331 + t277)) * m(3);
t8 = -t20 * t353 - t104 * t11 - t156 * t14 + 0.2e1 * ((t139 * t205 + t177) * t211 + t272 + t366) * t53 + (((t194 * t237 - t346) * rSges(3,1) + t229 * rSges(3,2)) * t211 + t205 * (rSges(3,1) * t229 - t237 * t331 + t277)) * m(3);
t7 = -t19 * t354 - t103 * t10 - t156 * t13 + 0.2e1 * ((t139 * t203 + t177) * t209 + t272 + t367) * t52 + (((t194 * t239 - t346) * rSges(3,1) + t230 * rSges(3,2)) * t209 + t203 * (rSges(3,1) * t230 - t239 * t331 + t277)) * m(3);
t6 = -t21 * t322 - t84 * t12 - t105 * t15 - 0.4e1 * ((t207 * t357 + t213 * t356) * t63 + ((t207 * t358 + t177) * t213 + t355 + t365) * t57) * t63 + (-t117 * t370 + t157 * t234) * t214 - (-t117 * t234 - t157 * t370) * t208;
t5 = -t20 * t323 - t83 * t11 - t104 * t14 - 0.4e1 * ((t205 * t357 + t211 * t356) * t62 + ((t205 * t358 + t177) * t211 + t355 + t366) * t56) * t62 + (-t116 * t369 + t157 * t236) * t212 - (-t116 * t236 - t157 * t369) * t206;
t4 = -t19 * t324 - t82 * t10 - t103 * t13 - 0.4e1 * ((t203 * t357 + t209 * t356) * t61 + ((t203 * t358 + t177) * t209 + t355 + t367) * t55) * t61 + (-t115 * t368 + t157 * t238) * t210 - (-t115 * t238 - t157 * t368) * t204;
t3 = -t12 * t322 - t348 + (-t196 * t325 - t328) * m(3) + t332;
t2 = -t11 * t323 - t348 + (-t196 * t326 - t329) * m(3) + t333;
t1 = -t10 * t324 - t348 + (-t196 * t327 - t330) * m(3) + t334;
t16 = [-m(4) * g(1) + (t66 * t3 + t72 * t6) * t87 + (t65 * t2 + t71 * t5) * t86 + (t64 * t1 + t70 * t4) * t85 + (t9 * t336 + t8 * t337 + t7 * t338) * t226 + ((t30 * t75 + t36 * t69) * t87 + (t29 * t74 + t35 * t68) * t86 + (t28 * t73 + t34 * t67) * t85 + (t339 * t42 + t340 * t41 + t341 * t40) * t226) * t201 + t264 * t200 + (m(4) + (t30 * t72 + t36 * t66) * t87 + (t29 * t71 + t35 * t65) * t86 + (t28 * t70 + t34 * t64) * t85 + (t336 * t42 + t337 * t41 + t338 * t40) * t226) * t202; -m(4) * g(2) + (t69 * t3 + t75 * t6) * t87 + (t68 * t2 + t74 * t5) * t86 + (t67 * t1 + t73 * t4) * t85 + (t9 * t339 + t8 * t340 + t7 * t341) * t226 + ((t33 * t72 + t39 * t66) * t87 + (t32 * t71 + t38 * t65) * t86 + (t31 * t70 + t37 * t64) * t85 + (t336 * t45 + t337 * t44 + t338 * t43) * t226) * t202 + t263 * t200 + (m(4) + (t33 * t75 + t39 * t69) * t87 + (t32 * t74 + t38 * t68) * t86 + (t31 * t73 + t37 * t67) * t85 + (t339 * t45 + t340 * t44 + t341 * t43) * t226) * t201; t264 * t202 + t263 * t201 + (-t97 * t10 - t98 * t11 - t99 * t12) * t194 + (-t330 - t329 - t328 + (-t325 - t326 - t327) * t196) * m(3) + t332 + t333 + t334 + (t200 - g(3)) * (0.3e1 * t189 + m(4));];
tauX  = t16;
