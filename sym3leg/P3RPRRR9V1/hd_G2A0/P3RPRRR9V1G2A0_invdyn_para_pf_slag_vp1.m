% Calculate vector of inverse dynamics forces for parallel robot
% P3RPRRR9V1G2A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-08-06 18:53
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRRR9V1G2A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:51:17
% EndTime: 2020-08-06 18:51:23
% DurationCPUTime: 6.40s
% Computational Cost: add. (23286->541), mult. (31119->893), div. (4275->10), fcn. (23826->40), ass. (0->365)
t247 = xDP(2);
t217 = cos(pkin(7));
t233 = cos(qJ(3,3));
t216 = sin(pkin(7));
t227 = sin(qJ(3,3));
t333 = t216 * t227;
t110 = 0.1e1 / (t217 * t233 - t333);
t394 = pkin(5) + qJ(2,3);
t202 = pkin(6) + t394;
t185 = 0.1e1 / t202;
t360 = t110 * t185;
t258 = -pkin(3) / 0.2e1;
t399 = t233 * pkin(2);
t212 = t233 ^ 2;
t402 = t212 * pkin(3);
t126 = t402 + t399 / 0.2e1 + t258;
t221 = legFrame(3,2);
t188 = sin(t221);
t191 = cos(t221);
t228 = sin(qJ(1,3));
t234 = cos(qJ(1,3));
t298 = pkin(1) * t228 - t202 * t234;
t305 = t228 * t333;
t422 = 0.2e1 * pkin(3);
t276 = pkin(2) * t305 + (t305 * t422 - t298) * t233;
t345 = t188 * t228;
t259 = pkin(2) / 0.2e1;
t351 = (pkin(3) * t233 + t259) * t227;
t184 = pkin(1) * t216;
t354 = (-pkin(3) * t227 + t184) * t233;
t420 = 0.2e1 * t217 ^ 2;
t76 = t298 * t333 + (t212 - 0.1e1) * t228 * pkin(3);
t91 = pkin(1) * t227 + (-pkin(3) + t399 + 0.2e1 * t402) * t216;
t49 = (-t126 * t345 + t191 * t351) * t420 + (t276 * t188 + t191 * t91) * t217 + t76 * t188 + t191 * t354;
t40 = t49 * t247 * t360;
t248 = xDP(1);
t339 = t191 * t228;
t50 = (t126 * t339 + t188 * t351) * t420 + (t188 * t91 - t276 * t191) * t217 - t76 * t191 + t188 * t354;
t41 = t50 * t248 * t360;
t246 = xDP(3);
t183 = t217 * pkin(2);
t161 = t183 + pkin(1);
t209 = pkin(7) + qJ(3,3);
t180 = cos(t209);
t405 = pkin(3) * t180;
t117 = t161 + t405;
t88 = t117 * t234 + t202 * t228;
t73 = t88 * t185 * t246;
t28 = t41 + t40 + t73;
t438 = 0.2e1 * t28;
t235 = cos(qJ(3,2));
t229 = sin(qJ(3,2));
t332 = t216 * t229;
t111 = 0.1e1 / (t217 * t235 - t332);
t395 = pkin(5) + qJ(2,2);
t203 = pkin(6) + t395;
t186 = 0.1e1 / t203;
t358 = t111 * t186;
t398 = t235 * pkin(2);
t213 = t235 ^ 2;
t401 = t213 * pkin(3);
t127 = t401 + t398 / 0.2e1 + t258;
t222 = legFrame(2,2);
t189 = sin(t222);
t192 = cos(t222);
t230 = sin(qJ(1,2));
t236 = cos(qJ(1,2));
t297 = pkin(1) * t230 - t203 * t236;
t304 = t230 * t332;
t275 = pkin(2) * t304 + (t304 * t422 - t297) * t235;
t343 = t189 * t230;
t350 = (pkin(3) * t235 + t259) * t229;
t353 = (-pkin(3) * t229 + t184) * t235;
t77 = t297 * t332 + (t213 - 0.1e1) * t230 * pkin(3);
t92 = pkin(1) * t229 + (-pkin(3) + t398 + 0.2e1 * t401) * t216;
t51 = (-t127 * t343 + t192 * t350) * t420 + (t275 * t189 + t192 * t92) * t217 + t77 * t189 + t192 * t353;
t42 = t51 * t247 * t358;
t337 = t192 * t230;
t52 = (t127 * t337 + t189 * t350) * t420 + (t189 * t92 - t275 * t192) * t217 - t77 * t192 + t189 * t353;
t43 = t52 * t248 * t358;
t210 = pkin(7) + qJ(3,2);
t181 = cos(t210);
t404 = pkin(3) * t181;
t118 = t161 + t404;
t89 = t118 * t236 + t203 * t230;
t74 = t89 * t186 * t246;
t29 = t43 + t42 + t74;
t437 = 0.2e1 * t29;
t237 = cos(qJ(3,1));
t231 = sin(qJ(3,1));
t331 = t216 * t231;
t112 = 0.1e1 / (t217 * t237 - t331);
t396 = pkin(5) + qJ(2,1);
t204 = pkin(6) + t396;
t187 = 0.1e1 / t204;
t356 = t112 * t187;
t397 = t237 * pkin(2);
t214 = t237 ^ 2;
t400 = t214 * pkin(3);
t128 = t400 + t397 / 0.2e1 + t258;
t223 = legFrame(1,2);
t190 = sin(t223);
t193 = cos(t223);
t232 = sin(qJ(1,1));
t238 = cos(qJ(1,1));
t296 = pkin(1) * t232 - t204 * t238;
t303 = t232 * t331;
t274 = pkin(2) * t303 + (t303 * t422 - t296) * t237;
t341 = t190 * t232;
t349 = (pkin(3) * t237 + t259) * t231;
t352 = (-pkin(3) * t231 + t184) * t237;
t78 = t296 * t331 + (t214 - 0.1e1) * t232 * pkin(3);
t93 = pkin(1) * t231 + (-pkin(3) + t397 + 0.2e1 * t400) * t216;
t53 = (-t128 * t341 + t193 * t349) * t420 + (t274 * t190 + t193 * t93) * t217 + t78 * t190 + t193 * t352;
t44 = t53 * t247 * t356;
t335 = t193 * t232;
t54 = (t128 * t335 + t190 * t349) * t420 + (t190 * t93 - t274 * t193) * t217 - t78 * t193 + t190 * t352;
t45 = t54 * t248 * t356;
t211 = pkin(7) + qJ(3,1);
t182 = cos(t211);
t403 = pkin(3) * t182;
t119 = t161 + t403;
t90 = t119 * t238 + t204 * t232;
t75 = t90 * t187 * t246;
t30 = t45 + t44 + t75;
t436 = 0.2e1 * t30;
t435 = t231 * rSges(3,1) + t237 * rSges(3,2);
t434 = t229 * rSges(3,1) + t235 * rSges(3,2);
t433 = t227 * rSges(3,1) + t233 * rSges(3,2);
t432 = 0.2e1 * pkin(1);
t168 = 0.1e1 / t180;
t174 = sin(t209);
t94 = t174 * t191 - t180 * t345;
t95 = t174 * t188 + t180 * t339;
t58 = (t234 * t246 + (t247 * t94 + t248 * t95) * t168) * t185;
t431 = t58 / 0.2e1;
t169 = 0.1e1 / t181;
t175 = sin(t210);
t96 = t175 * t192 - t181 * t343;
t97 = t175 * t189 + t181 * t337;
t59 = (t236 * t246 + (t247 * t96 + t248 * t97) * t169) * t186;
t430 = t59 / 0.2e1;
t170 = 0.1e1 / t182;
t176 = sin(t211);
t98 = t176 * t193 - t182 * t341;
t99 = t176 * t190 + t182 * t335;
t60 = (t238 * t246 + (t247 * t98 + t248 * t99) * t170) * t187;
t429 = t60 / 0.2e1;
t428 = -0.2e1 * t183;
t427 = t432 / 0.2e1;
t323 = 2 * m(3);
t426 = (rSges(3,1) * t176 + rSges(3,2) * t182) * t323;
t425 = (rSges(3,1) * t175 + rSges(3,2) * t181) * t323;
t424 = (rSges(3,1) * t174 + rSges(3,2) * t180) * t323;
t243 = rSges(2,2) * m(2);
t415 = m(2) * rSges(2,1);
t416 = pkin(2) * m(3);
t272 = -t216 * t243 + (t415 + t416) * t217;
t421 = 4 * rSges(2,3);
t419 = -4 * pkin(5) - 4 * pkin(6);
t418 = m(2) / 0.2e1;
t417 = m(3) / 0.2e1;
t414 = m(3) * rSges(3,1);
t413 = m(3) * rSges(3,2);
t263 = rSges(3,2) ^ 2;
t265 = rSges(3,1) ^ 2;
t138 = m(3) * (-t263 + t265) - Icges(3,1) + Icges(3,2);
t412 = t138 / 0.2e1;
t411 = m(2) * (rSges(2,3) + qJ(2,3));
t410 = m(2) * (rSges(2,3) + qJ(2,2));
t409 = m(2) * (rSges(2,3) + qJ(2,1));
t198 = (rSges(3,3) + t394);
t408 = m(3) * t198;
t199 = (rSges(3,3) + t395);
t407 = m(3) * t199;
t200 = (rSges(3,3) + t396);
t406 = m(3) * t200;
t393 = t110 * t49;
t392 = t110 * t50;
t391 = t111 * t51;
t390 = t111 * t52;
t389 = t112 * t53;
t388 = t112 * t54;
t162 = 0.2e1 * t209;
t147 = sin(t162);
t387 = t147 * t58;
t163 = 0.2e1 * t210;
t148 = sin(t163);
t386 = t148 * t59;
t164 = 0.2e1 * t211;
t149 = sin(t164);
t385 = t149 * t60;
t268 = 0.1e1 / pkin(3);
t85 = (t188 * t248 + t191 * t247) * t268 * t168;
t384 = t85 ^ 2 * t168;
t383 = t168 * t94;
t382 = t168 * t95;
t86 = (t189 * t248 + t192 * t247) * t268 * t169;
t381 = t86 ^ 2 * t169;
t380 = t169 * t96;
t379 = t169 * t97;
t87 = (t190 * t248 + t193 * t247) * t268 * t170;
t378 = t87 ^ 2 * t170;
t377 = t170 * t98;
t376 = t170 * t99;
t375 = t174 * t85;
t374 = t175 * t86;
t373 = t176 * t87;
t372 = t185 * t94;
t371 = t185 * t95;
t370 = t186 * t96;
t369 = t186 * t97;
t368 = t187 * t98;
t367 = t187 * t99;
t250 = m(2) + m(3);
t359 = t110 * t250;
t357 = t111 * t250;
t355 = t112 * t250;
t150 = cos(t162);
t166 = rSges(3,1) * t413 - Icges(3,4);
t348 = t166 * t150;
t151 = cos(t163);
t347 = t166 * t151;
t152 = cos(t164);
t346 = t166 * t152;
t344 = t188 * t268;
t342 = t189 * t268;
t340 = t190 * t268;
t338 = t191 * t268;
t336 = t192 * t268;
t334 = t193 * t268;
t257 = 0.2e1 * pkin(7);
t206 = t257 + qJ(3,3);
t171 = sin(t206);
t330 = t171 + t227;
t207 = t257 + qJ(3,2);
t172 = sin(t207);
t329 = t172 + t229;
t208 = t257 + qJ(3,1);
t173 = sin(t208);
t328 = t173 + t231;
t177 = cos(t206);
t327 = t177 + t233;
t178 = cos(t207);
t326 = t178 + t235;
t179 = cos(t208);
t325 = t179 + t237;
t324 = t263 + t265;
t279 = -t174 * t413 + t180 * t414;
t289 = pkin(1) * t250 + t272;
t70 = t279 + t289;
t316 = t70 * t360;
t278 = -t175 * t413 + t181 * t414;
t71 = t278 + t289;
t315 = t71 * t358;
t277 = -t176 * t413 + t182 * t414;
t72 = t277 + t289;
t314 = t72 * t356;
t313 = t174 * t384;
t312 = t175 * t381;
t311 = t176 * t378;
t140 = -rSges(3,2) * t408 + Icges(3,6);
t292 = rSges(3,1) * t408 - Icges(3,5);
t79 = t140 * t180 - t292 * t174;
t310 = t234 * t268 * t79;
t141 = -rSges(3,2) * t407 + Icges(3,6);
t291 = rSges(3,1) * t407 - Icges(3,5);
t80 = t141 * t181 - t291 * t175;
t309 = t236 * t268 * t80;
t142 = -rSges(3,2) * t406 + Icges(3,6);
t290 = rSges(3,1) * t406 - Icges(3,5);
t81 = t142 * t182 - t290 * t176;
t308 = t238 * t268 * t81;
t307 = t414 * t432;
t306 = -0.2e1 * pkin(1) * t413;
t165 = m(1) * rSges(1,1) + m(2) * pkin(1);
t302 = -m(3) * pkin(1) - t165;
t270 = pkin(1) ^ 2;
t255 = 0.2e1 * t270;
t264 = rSges(2,2) ^ 2;
t266 = rSges(2,1) ^ 2;
t295 = (2 * rSges(2,3) ^ 2) + t255 + t264 + t266;
t269 = pkin(2) ^ 2;
t294 = t255 + t269 + t324;
t123 = g(1) * t191 - g(2) * t188;
t285 = g(3) * t234 + t123 * t228;
t124 = g(1) * t192 - g(2) * t189;
t284 = g(3) * t236 + t124 * t230;
t125 = g(1) * t193 - g(2) * t190;
t283 = g(3) * t238 + t125 * t232;
t282 = -(-t415 + (-rSges(3,1) * t233 + rSges(3,2) * t227 - pkin(2)) * m(3)) * t217 - (t433 * m(3) + t243) * t216;
t281 = -(-t415 + (-rSges(3,1) * t235 + rSges(3,2) * t229 - pkin(2)) * m(3)) * t217 - (t434 * m(3) + t243) * t216;
t280 = -(-t415 + (-rSges(3,1) * t237 + rSges(3,2) * t231 - pkin(2)) * m(3)) * t217 - (t435 * m(3) + t243) * t216;
t201 = cos(t257);
t267 = pkin(3) ^ 2;
t273 = -(2 * pkin(6) ^ 2) - t269 * t201 - t267 - t269 - 0.2e1 * t270 + ((-4 * pkin(6) - 2 * pkin(5)) * pkin(5));
t271 = Icges(1,3) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + (-rSges(2,1) * t243 + Icges(2,4)) * sin(t257) + Icges(2,1) / 0.2e1 + Icges(3,1) / 0.2e1 + Icges(2,2) / 0.2e1 + Icges(3,2) / 0.2e1 + (m(3) * t269 + (-t264 + t266) * m(2) - Icges(2,1) + Icges(2,2)) * t201 / 0.2e1;
t262 = qJ(2,1) ^ 2;
t261 = qJ(2,2) ^ 2;
t260 = qJ(2,3) ^ 2;
t252 = pkin(1) * g(3);
t244 = m(1) * rSges(1,2);
t226 = xDDP(1);
t225 = xDDP(2);
t224 = xDDP(3);
t156 = t165 * g(3);
t146 = t324 * m(3) + Icges(3,3);
t145 = t244 - t409;
t144 = t244 - t410;
t143 = t244 - t411;
t132 = t406 + t409;
t131 = t407 + t410;
t130 = t408 + t411;
t122 = g(1) * t190 + g(2) * t193;
t121 = g(1) * t189 + g(2) * t192;
t120 = g(1) * t188 + g(2) * t191;
t69 = (-t238 * t72 + t250 * t90) * t187;
t68 = (-t236 * t71 + t250 * t89) * t186;
t67 = (-t234 * t70 + t250 * t88) * t185;
t66 = (t146 * t340 + t81 * t367) * t170;
t65 = (t146 * t334 + t81 * t368) * t170;
t64 = (t146 * t342 + t80 * t369) * t169;
t63 = (t146 * t336 + t80 * t370) * t169;
t62 = (t146 * t344 + t79 * t371) * t168;
t61 = (t146 * t338 + t79 * t372) * t168;
t57 = pkin(1) * t60;
t56 = pkin(1) * t59;
t55 = pkin(1) * t58;
t48 = t152 * t412 + ((2 * t200 ^ 2) + t294) * t417 + ((qJ(2,1) * t421) + (2 * t262) + t295) * t418 - t166 * t149 + (t325 * rSges(3,1) - t328 * rSges(3,2)) * t416 + (t272 + t277) * t432 + t271;
t47 = t151 * t412 + ((2 * t199 ^ 2) + t294) * t417 + ((qJ(2,2) * t421) + (2 * t261) + t295) * t418 - t166 * t148 + (t326 * rSges(3,1) - t329 * rSges(3,2)) * t416 + (t272 + t278) * t432 + t271;
t46 = t150 * t412 + ((2 * t198 ^ 2) + t294) * t417 + ((qJ(2,3) * t421) + (2 * t260) + t295) * t418 - t166 * t147 + (t327 * rSges(3,1) - t330 * rSges(3,2)) * t416 + (t272 + t279) * t432 + t271;
t39 = (t238 * t48 - t72 * t90) * t187;
t38 = (t236 * t47 - t71 * t89) * t186;
t37 = (t234 * t46 - t70 * t88) * t185;
t36 = (t54 * t355 - t72 * t376) * t187;
t35 = (t53 * t355 - t72 * t377) * t187;
t34 = (t52 * t357 - t71 * t379) * t186;
t33 = (t51 * t357 - t71 * t380) * t186;
t32 = (t50 * t359 - t70 * t382) * t185;
t31 = (t49 * t359 - t70 * t383) * t185;
t27 = -t54 * t314 + (t81 * t340 + t48 * t367) * t170;
t26 = -t53 * t314 + (t81 * t334 + t48 * t368) * t170;
t25 = -t52 * t315 + (t80 * t342 + t47 * t369) * t169;
t24 = -t51 * t315 + (t80 * t336 + t47 * t370) * t169;
t23 = -t50 * t316 + (t79 * t344 + t46 * t371) * t168;
t22 = -t49 * t316 + (t79 * t338 + t46 * t372) * t168;
t18 = t57 - t45 / 0.2e1 - t44 / 0.2e1 - t75 / 0.2e1;
t17 = t56 - t43 / 0.2e1 - t42 / 0.2e1 - t74 / 0.2e1;
t16 = t55 - t41 / 0.2e1 - t40 / 0.2e1 - t73 / 0.2e1;
t15 = (-pkin(3) * t378 + (-t57 + t436 + (-t183 - t403) * t60) * t60) * t187;
t14 = (-pkin(3) * t381 + (-t56 + t437 + (-t183 - t404) * t59) * t59) * t186;
t13 = (-pkin(3) * t384 + (-t55 + t438 + (-t183 - t405) * t58) * t58) * t185;
t12 = ((t18 * t428 + ((qJ(2,1) * t419) - t267 * t152 - (2 * t262) + t273) * t429 + (t119 + t427) * t30) * t60 + ((-0.2e1 * t18 * t182 + t204 * t373 + (-pkin(2) * t179 - t397) * t60) * t60 - (-t204 * t385 / 0.2e1 + t87 * t119) * t170 * t87) * pkin(3)) * t187;
t11 = ((t17 * t428 + ((qJ(2,2) * t419) - t267 * t151 - (2 * t261) + t273) * t430 + (t118 + t427) * t29) * t59 + ((-0.2e1 * t17 * t181 + t203 * t374 + (-pkin(2) * t178 - t398) * t59) * t59 - (-t203 * t386 / 0.2e1 + t86 * t118) * t169 * t86) * pkin(3)) * t186;
t10 = ((t16 * t428 + ((qJ(2,3) * t419) - t267 * t150 - (2 * t260) + t273) * t431 + (t117 + t427) * t28) * t58 + ((-0.2e1 * t16 * t180 + t202 * t375 + (-pkin(2) * t177 - t399) * t58) * t58 - (-t202 * t387 / 0.2e1 + t85 * t117) * t168 * t85) * pkin(3)) * t185;
t9 = -t81 * t15 + t146 * t311 + ((t57 - 0.2e1 * t45 - 0.2e1 * t44 - 0.2e1 * t75) * t426 + (t138 * t149 + 0.2e1 * t346 + (t328 * rSges(3,1) + t325 * rSges(3,2)) * t416) * t60) * t429 + m(3) * ((-rSges(3,1) * t122 + t283 * rSges(3,2)) * t182 + t176 * (t283 * rSges(3,1) + rSges(3,2) * t122));
t8 = -t80 * t14 + t146 * t312 + ((t56 - 0.2e1 * t43 - 0.2e1 * t42 - 0.2e1 * t74) * t425 + (t138 * t148 + 0.2e1 * t347 + (t329 * rSges(3,1) + t326 * rSges(3,2)) * t416) * t59) * t430 + m(3) * ((-rSges(3,1) * t121 + t284 * rSges(3,2)) * t181 + t175 * (t284 * rSges(3,1) + rSges(3,2) * t121));
t7 = -t79 * t13 + t146 * t313 + ((t55 - 0.2e1 * t41 - 0.2e1 * t40 - 0.2e1 * t73) * t424 + (t138 * t147 + 0.2e1 * t348 + (t330 * rSges(3,1) + t327 * rSges(3,2)) * t416) * t58) * t431 + m(3) * ((-rSges(3,1) * t120 + t285 * rSges(3,2)) * t180 + t174 * (t285 * rSges(3,1) + rSges(3,2) * t120));
t6 = t72 * t15 - t60 * (t132 * t60 - t87 * t426) + (-g(3) * t232 + t125 * t238 - t12) * t250;
t5 = t71 * t14 - t59 * (t131 * t59 - t86 * t425) + (-g(3) * t230 + t124 * t236 - t11) * t250;
t4 = t70 * t13 - t58 * (t130 * t58 - t85 * t424) + (-g(3) * t228 + t123 * t234 - t10) * t250;
t3 = -t48 * t15 + t72 * t12 + t81 * t311 + ((t145 - t406) * g(3) + (-t280 + t302) * t125) * t238 + ((-t125 * t200 + t252) * m(3) + t156 + t145 * t125 + t280 * g(3)) * t232 + (-t290 * t182 * t87 - t138 * t385 - t142 * t373) * t87 + (t132 * t436 - t307 * t373 + (t306 * t182 - 0.2e1 * t346 + (-rSges(3,1) * t173 - rSges(3,2) * t179 - t435) * t416) * t87) * t60;
t2 = -t47 * t14 + t71 * t11 + t80 * t312 + ((t144 - t407) * g(3) + (-t281 + t302) * t124) * t236 + ((-t124 * t199 + t252) * m(3) + t156 + t144 * t124 + t281 * g(3)) * t230 + (-t291 * t181 * t86 - t138 * t386 - t141 * t374) * t86 + (t131 * t437 - t307 * t374 + (t306 * t181 - 0.2e1 * t347 + (-rSges(3,1) * t172 - rSges(3,2) * t178 - t434) * t416) * t86) * t59;
t1 = -t46 * t13 + t70 * t10 + t79 * t313 + ((t143 - t408) * g(3) + (-t282 + t302) * t123) * t234 + ((-t123 * t198 + t252) * m(3) + t156 + t143 * t123 + t282 * g(3)) * t228 + (-t292 * t180 * t85 - t138 * t387 - t140 * t375) * t85 + (t130 * t438 - t307 * t375 + (t306 * t180 - 0.2e1 * t348 + (-rSges(3,1) * t171 - rSges(3,2) * t177 - t433) * t416) * t85) * t58;
t19 = [(-g(1) + t226) * m(4) + ((t27 * t376 + t36 * t388) * t226 + (t27 * t377 + t36 * t389) * t225 + (t238 * t27 + t36 * t90) * t224 + t3 * t376 + t6 * t388) * t187 + ((t193 * t225 * t66 + (t226 * t66 + t9) * t190) * t170 + (t192 * t225 * t64 + (t226 * t64 + t8) * t189) * t169 + (t191 * t225 * t62 + (t226 * t62 + t7) * t188) * t168) * t268 + ((t25 * t379 + t34 * t390) * t226 + (t25 * t380 + t34 * t391) * t225 + (t236 * t25 + t34 * t89) * t224 + t2 * t379 + t5 * t390) * t186 + ((t23 * t382 + t32 * t392) * t226 + (t23 * t383 + t32 * t393) * t225 + (t23 * t234 + t32 * t88) * t224 + t1 * t382 + t4 * t392) * t185; (-g(2) + t225) * m(4) + ((t26 * t376 + t35 * t388) * t226 + (t26 * t377 + t35 * t389) * t225 + (t238 * t26 + t35 * t90) * t224 + t3 * t377 + t6 * t389) * t187 + ((t190 * t226 * t65 + (t225 * t65 + t9) * t193) * t170 + (t189 * t226 * t63 + (t225 * t63 + t8) * t192) * t169 + (t188 * t226 * t61 + (t225 * t61 + t7) * t191) * t168) * t268 + ((t24 * t379 + t33 * t390) * t226 + (t24 * t380 + t33 * t391) * t225 + (t236 * t24 + t33 * t89) * t224 + t2 * t380 + t5 * t391) * t186 + ((t22 * t382 + t31 * t392) * t226 + (t22 * t383 + t31 * t393) * t225 + (t22 * t234 + t31 * t88) * t224 + t1 * t383 + t4 * t393) * t185; (-g(3) + t224) * m(4) + ((t238 * t39 + t69 * t90) * t224 + t238 * t3 + t90 * t6 + (t225 * t53 + t226 * t54) * t69 * t112 + ((t190 * t308 + t39 * t99) * t226 + (t193 * t308 + t39 * t98) * t225) * t170) * t187 + ((t236 * t38 + t68 * t89) * t224 + t236 * t2 + t89 * t5 + (t225 * t51 + t226 * t52) * t68 * t111 + ((t189 * t309 + t38 * t97) * t226 + (t192 * t309 + t38 * t96) * t225) * t169) * t186 + ((t234 * t37 + t67 * t88) * t224 + t234 * t1 + t88 * t4 + (t225 * t49 + t226 * t50) * t67 * t110 + ((t188 * t310 + t37 * t95) * t226 + (t191 * t310 + t37 * t94) * t225) * t168) * t185;];
tauX  = t19;
