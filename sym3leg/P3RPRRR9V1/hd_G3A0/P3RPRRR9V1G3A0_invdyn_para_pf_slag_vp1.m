% Calculate vector of inverse dynamics forces for parallel robot
% P3RPRRR9V1G3A0
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
% Datum: 2020-08-06 18:58
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRRR9V1G3A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G3A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR9V1G3A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR9V1G3A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G3A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR9V1G3A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G3A0_invdyn_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR9V1G3A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR9V1G3A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRRR9V1G3A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G3A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G3A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:56:15
% EndTime: 2020-08-06 18:56:22
% DurationCPUTime: 6.42s
% Computational Cost: add. (23286->541), mult. (31119->893), div. (4275->10), fcn. (23826->40), ass. (0->365)
t250 = xDP(2);
t220 = cos(pkin(7));
t236 = cos(qJ(3,3));
t219 = sin(pkin(7));
t230 = sin(qJ(3,3));
t336 = t219 * t230;
t110 = 0.1e1 / (t220 * t236 - t336);
t397 = pkin(5) + qJ(2,3);
t205 = pkin(6) + t397;
t188 = 0.1e1 / t205;
t363 = t110 * t188;
t261 = -pkin(3) / 0.2e1;
t402 = t236 * pkin(2);
t215 = t236 ^ 2;
t405 = t215 * pkin(3);
t126 = t405 + t402 / 0.2e1 + t261;
t224 = legFrame(3,2);
t191 = sin(t224);
t194 = cos(t224);
t237 = cos(qJ(1,3));
t305 = t237 * t336;
t231 = sin(qJ(1,3));
t327 = pkin(1) * t237 + t231 * t205;
t425 = 0.2e1 * pkin(3);
t279 = pkin(2) * t305 + (t305 * t425 - t327) * t236;
t348 = t191 * t237;
t262 = pkin(2) / 0.2e1;
t354 = (pkin(3) * t236 + t262) * t230;
t187 = pkin(1) * t219;
t357 = (-pkin(3) * t230 + t187) * t236;
t423 = 0.2e1 * t220 ^ 2;
t76 = t327 * t336 + (t215 - 0.1e1) * t237 * pkin(3);
t91 = pkin(1) * t230 + (-pkin(3) + t402 + 0.2e1 * t405) * t219;
t49 = (-t126 * t348 + t194 * t354) * t423 + (t279 * t191 + t194 * t91) * t220 + t76 * t191 + t194 * t357;
t40 = t49 * t250 * t363;
t251 = xDP(1);
t342 = t194 * t237;
t50 = (t126 * t342 + t191 * t354) * t423 + (t191 * t91 - t279 * t194) * t220 - t76 * t194 + t191 * t357;
t41 = t50 * t251 * t363;
t249 = xDP(3);
t186 = t220 * pkin(2);
t164 = t186 + pkin(1);
t212 = pkin(7) + qJ(3,3);
t183 = cos(t212);
t408 = pkin(3) * t183;
t117 = t164 + t408;
t88 = -t117 * t231 + t205 * t237;
t73 = t88 * t188 * t249;
t28 = t41 + t40 + t73;
t441 = 0.2e1 * t28;
t238 = cos(qJ(3,2));
t232 = sin(qJ(3,2));
t335 = t219 * t232;
t111 = 0.1e1 / (t220 * t238 - t335);
t398 = pkin(5) + qJ(2,2);
t206 = pkin(6) + t398;
t189 = 0.1e1 / t206;
t361 = t111 * t189;
t401 = t238 * pkin(2);
t216 = t238 ^ 2;
t404 = t216 * pkin(3);
t127 = t404 + t401 / 0.2e1 + t261;
t225 = legFrame(2,2);
t192 = sin(t225);
t195 = cos(t225);
t239 = cos(qJ(1,2));
t304 = t239 * t335;
t233 = sin(qJ(1,2));
t326 = pkin(1) * t239 + t233 * t206;
t278 = pkin(2) * t304 + (t304 * t425 - t326) * t238;
t346 = t192 * t239;
t353 = (pkin(3) * t238 + t262) * t232;
t356 = (-pkin(3) * t232 + t187) * t238;
t77 = t326 * t335 + (t216 - 0.1e1) * t239 * pkin(3);
t92 = pkin(1) * t232 + (-pkin(3) + t401 + 0.2e1 * t404) * t219;
t51 = (-t127 * t346 + t195 * t353) * t423 + (t278 * t192 + t195 * t92) * t220 + t77 * t192 + t195 * t356;
t42 = t51 * t250 * t361;
t340 = t195 * t239;
t52 = (t127 * t340 + t192 * t353) * t423 + (t192 * t92 - t278 * t195) * t220 - t77 * t195 + t192 * t356;
t43 = t52 * t251 * t361;
t213 = pkin(7) + qJ(3,2);
t184 = cos(t213);
t407 = pkin(3) * t184;
t118 = t164 + t407;
t89 = -t118 * t233 + t206 * t239;
t74 = t89 * t189 * t249;
t29 = t43 + t42 + t74;
t440 = 0.2e1 * t29;
t240 = cos(qJ(3,1));
t234 = sin(qJ(3,1));
t334 = t219 * t234;
t112 = 0.1e1 / (t220 * t240 - t334);
t399 = pkin(5) + qJ(2,1);
t207 = pkin(6) + t399;
t190 = 0.1e1 / t207;
t359 = t112 * t190;
t400 = t240 * pkin(2);
t217 = t240 ^ 2;
t403 = t217 * pkin(3);
t128 = t403 + t400 / 0.2e1 + t261;
t226 = legFrame(1,2);
t193 = sin(t226);
t196 = cos(t226);
t241 = cos(qJ(1,1));
t303 = t241 * t334;
t235 = sin(qJ(1,1));
t325 = pkin(1) * t241 + t235 * t207;
t277 = pkin(2) * t303 + (t303 * t425 - t325) * t240;
t344 = t193 * t241;
t352 = (pkin(3) * t240 + t262) * t234;
t355 = (-pkin(3) * t234 + t187) * t240;
t78 = t325 * t334 + (t217 - 0.1e1) * t241 * pkin(3);
t93 = pkin(1) * t234 + (-pkin(3) + t400 + 0.2e1 * t403) * t219;
t53 = (-t128 * t344 + t196 * t352) * t423 + (t277 * t193 + t196 * t93) * t220 + t78 * t193 + t196 * t355;
t44 = t53 * t250 * t359;
t338 = t196 * t241;
t54 = (t128 * t338 + t193 * t352) * t423 + (t193 * t93 - t277 * t196) * t220 - t78 * t196 + t193 * t355;
t45 = t54 * t251 * t359;
t214 = pkin(7) + qJ(3,1);
t185 = cos(t214);
t406 = pkin(3) * t185;
t119 = t164 + t406;
t90 = -t119 * t235 + t207 * t241;
t75 = t90 * t190 * t249;
t30 = t45 + t44 + t75;
t439 = 0.2e1 * t30;
t438 = t234 * rSges(3,1) + t240 * rSges(3,2);
t437 = t232 * rSges(3,1) + t238 * rSges(3,2);
t436 = t230 * rSges(3,1) + t236 * rSges(3,2);
t435 = 0.2e1 * pkin(1);
t171 = 0.1e1 / t183;
t177 = sin(t212);
t94 = t177 * t194 - t183 * t348;
t95 = t177 * t191 + t183 * t342;
t58 = (-t231 * t249 + (t250 * t94 + t251 * t95) * t171) * t188;
t434 = t58 / 0.2e1;
t172 = 0.1e1 / t184;
t178 = sin(t213);
t96 = t178 * t195 - t184 * t346;
t97 = t178 * t192 + t184 * t340;
t59 = (-t233 * t249 + (t250 * t96 + t251 * t97) * t172) * t189;
t433 = t59 / 0.2e1;
t173 = 0.1e1 / t185;
t179 = sin(t214);
t98 = t179 * t196 - t185 * t344;
t99 = t179 * t193 + t185 * t338;
t60 = (-t235 * t249 + (t250 * t98 + t251 * t99) * t173) * t190;
t432 = t60 / 0.2e1;
t431 = -0.2e1 * t186;
t430 = t435 / 0.2e1;
t323 = 2 * m(3);
t429 = (rSges(3,1) * t179 + rSges(3,2) * t185) * t323;
t428 = (rSges(3,1) * t178 + rSges(3,2) * t184) * t323;
t427 = (rSges(3,1) * t177 + rSges(3,2) * t183) * t323;
t246 = rSges(2,2) * m(2);
t418 = m(2) * rSges(2,1);
t419 = pkin(2) * m(3);
t275 = -t219 * t246 + (t418 + t419) * t220;
t424 = 4 * rSges(2,3);
t422 = -4 * pkin(5) - 4 * pkin(6);
t421 = m(2) / 0.2e1;
t420 = m(3) / 0.2e1;
t417 = m(3) * rSges(3,1);
t416 = m(3) * rSges(3,2);
t266 = rSges(3,2) ^ 2;
t268 = rSges(3,1) ^ 2;
t138 = m(3) * (-t266 + t268) - Icges(3,1) + Icges(3,2);
t415 = t138 / 0.2e1;
t414 = m(2) * (rSges(2,3) + qJ(2,3));
t413 = m(2) * (rSges(2,3) + qJ(2,2));
t412 = m(2) * (rSges(2,3) + qJ(2,1));
t201 = (rSges(3,3) + t397);
t411 = m(3) * t201;
t202 = (rSges(3,3) + t398);
t410 = m(3) * t202;
t203 = (rSges(3,3) + t399);
t409 = m(3) * t203;
t396 = t110 * t49;
t395 = t110 * t50;
t394 = t111 * t51;
t393 = t111 * t52;
t392 = t112 * t53;
t391 = t112 * t54;
t165 = 0.2e1 * t212;
t147 = sin(t165);
t390 = t147 * t58;
t166 = 0.2e1 * t213;
t148 = sin(t166);
t389 = t148 * t59;
t167 = 0.2e1 * t214;
t149 = sin(t167);
t388 = t149 * t60;
t271 = 0.1e1 / pkin(3);
t85 = (t191 * t251 + t194 * t250) * t271 * t171;
t387 = t85 ^ 2 * t171;
t386 = t171 * t94;
t385 = t171 * t95;
t86 = (t192 * t251 + t195 * t250) * t271 * t172;
t384 = t86 ^ 2 * t172;
t383 = t172 * t96;
t382 = t172 * t97;
t87 = (t193 * t251 + t196 * t250) * t271 * t173;
t381 = t87 ^ 2 * t173;
t380 = t173 * t98;
t379 = t173 * t99;
t378 = t177 * t85;
t377 = t178 * t86;
t376 = t179 * t87;
t375 = t188 * t94;
t374 = t188 * t95;
t373 = t189 * t96;
t372 = t189 * t97;
t371 = t190 * t98;
t370 = t190 * t99;
t253 = m(2) + m(3);
t362 = t110 * t253;
t360 = t111 * t253;
t358 = t112 * t253;
t150 = cos(t165);
t169 = rSges(3,1) * t416 - Icges(3,4);
t351 = t169 * t150;
t151 = cos(t166);
t350 = t169 * t151;
t152 = cos(t167);
t349 = t169 * t152;
t347 = t191 * t271;
t345 = t192 * t271;
t343 = t193 * t271;
t341 = t194 * t271;
t339 = t195 * t271;
t337 = t196 * t271;
t260 = 0.2e1 * pkin(7);
t209 = t260 + qJ(3,3);
t174 = sin(t209);
t333 = t174 + t230;
t210 = t260 + qJ(3,2);
t175 = sin(t210);
t332 = t175 + t232;
t211 = t260 + qJ(3,1);
t176 = sin(t211);
t331 = t176 + t234;
t180 = cos(t209);
t330 = t180 + t236;
t181 = cos(t210);
t329 = t181 + t238;
t182 = cos(t211);
t328 = t182 + t240;
t324 = t266 + t268;
t282 = -t177 * t416 + t183 * t417;
t292 = pkin(1) * t253 + t275;
t70 = t282 + t292;
t316 = t70 * t363;
t281 = -t178 * t416 + t184 * t417;
t71 = t281 + t292;
t315 = t71 * t361;
t280 = -t179 * t416 + t185 * t417;
t72 = t280 + t292;
t314 = t72 * t359;
t313 = t177 * t387;
t312 = t178 * t384;
t311 = t179 * t381;
t140 = -rSges(3,2) * t411 + Icges(3,6);
t295 = rSges(3,1) * t411 - Icges(3,5);
t79 = t140 * t183 - t295 * t177;
t310 = t231 * t271 * t79;
t141 = -rSges(3,2) * t410 + Icges(3,6);
t294 = rSges(3,1) * t410 - Icges(3,5);
t80 = t141 * t184 - t294 * t178;
t309 = t233 * t271 * t80;
t142 = -rSges(3,2) * t409 + Icges(3,6);
t293 = rSges(3,1) * t409 - Icges(3,5);
t81 = t142 * t185 - t293 * t179;
t308 = t235 * t271 * t81;
t307 = t417 * t435;
t306 = -0.2e1 * pkin(1) * t416;
t168 = m(1) * rSges(1,1) + m(2) * pkin(1);
t302 = m(3) * pkin(1) + t168;
t273 = pkin(1) ^ 2;
t258 = 0.2e1 * t273;
t267 = rSges(2,2) ^ 2;
t269 = rSges(2,1) ^ 2;
t298 = (2 * rSges(2,3) ^ 2) + t258 + t267 + t269;
t272 = pkin(2) ^ 2;
t297 = t258 + t272 + t324;
t123 = g(1) * t194 - g(2) * t191;
t288 = -g(3) * t231 + t123 * t237;
t124 = g(1) * t195 - g(2) * t192;
t287 = -g(3) * t233 + t124 * t239;
t125 = g(1) * t196 - g(2) * t193;
t286 = -g(3) * t235 + t125 * t241;
t285 = -(-t418 + (-rSges(3,1) * t236 + rSges(3,2) * t230 - pkin(2)) * m(3)) * t220 - (t436 * m(3) + t246) * t219;
t284 = -(-t418 + (-rSges(3,1) * t238 + rSges(3,2) * t232 - pkin(2)) * m(3)) * t220 - (t437 * m(3) + t246) * t219;
t283 = -(-t418 + (-rSges(3,1) * t240 + rSges(3,2) * t234 - pkin(2)) * m(3)) * t220 - (t438 * m(3) + t246) * t219;
t204 = cos(t260);
t270 = pkin(3) ^ 2;
t276 = -(2 * pkin(6) ^ 2) - t272 * t204 - t270 - t272 - 0.2e1 * t273 + ((-4 * pkin(6) - 2 * pkin(5)) * pkin(5));
t274 = Icges(1,3) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + (-rSges(2,1) * t246 + Icges(2,4)) * sin(t260) + Icges(2,1) / 0.2e1 + Icges(3,1) / 0.2e1 + Icges(2,2) / 0.2e1 + Icges(3,2) / 0.2e1 + (m(3) * t272 + (-t267 + t269) * m(2) - Icges(2,1) + Icges(2,2)) * t204 / 0.2e1;
t265 = qJ(2,1) ^ 2;
t264 = qJ(2,2) ^ 2;
t263 = qJ(2,3) ^ 2;
t255 = g(3) * pkin(1);
t247 = m(1) * rSges(1,2);
t229 = xDDP(1);
t228 = xDDP(2);
t227 = xDDP(3);
t156 = g(3) * t168;
t146 = t324 * m(3) + Icges(3,3);
t145 = t247 - t412;
t144 = t247 - t413;
t143 = t247 - t414;
t132 = t409 + t412;
t131 = t410 + t413;
t130 = t411 + t414;
t122 = g(1) * t193 + g(2) * t196;
t121 = g(1) * t192 + g(2) * t195;
t120 = g(1) * t191 + g(2) * t194;
t69 = (t235 * t72 + t253 * t90) * t190;
t68 = (t233 * t71 + t253 * t89) * t189;
t67 = (t231 * t70 + t253 * t88) * t188;
t66 = (t146 * t343 + t81 * t370) * t173;
t65 = (t146 * t337 + t81 * t371) * t173;
t64 = (t146 * t345 + t80 * t372) * t172;
t63 = (t146 * t339 + t80 * t373) * t172;
t62 = (t146 * t347 + t79 * t374) * t171;
t61 = (t146 * t341 + t79 * t375) * t171;
t57 = pkin(1) * t60;
t56 = pkin(1) * t59;
t55 = pkin(1) * t58;
t48 = t152 * t415 + ((2 * t203 ^ 2) + t297) * t420 + ((qJ(2,1) * t424) + (2 * t265) + t298) * t421 - t169 * t149 + (t328 * rSges(3,1) - t331 * rSges(3,2)) * t419 + (t275 + t280) * t435 + t274;
t47 = t151 * t415 + ((2 * t202 ^ 2) + t297) * t420 + ((qJ(2,2) * t424) + (2 * t264) + t298) * t421 - t169 * t148 + (t329 * rSges(3,1) - t332 * rSges(3,2)) * t419 + (t275 + t281) * t435 + t274;
t46 = t150 * t415 + ((2 * t201 ^ 2) + t297) * t420 + ((qJ(2,3) * t424) + (2 * t263) + t298) * t421 - t169 * t147 + (t330 * rSges(3,1) - t333 * rSges(3,2)) * t419 + (t275 + t282) * t435 + t274;
t39 = (-t235 * t48 - t72 * t90) * t190;
t38 = (-t233 * t47 - t71 * t89) * t189;
t37 = (-t231 * t46 - t70 * t88) * t188;
t36 = (t54 * t358 - t72 * t379) * t190;
t35 = (t53 * t358 - t72 * t380) * t190;
t34 = (t52 * t360 - t71 * t382) * t189;
t33 = (t51 * t360 - t71 * t383) * t189;
t32 = (t50 * t362 - t70 * t385) * t188;
t31 = (t49 * t362 - t70 * t386) * t188;
t27 = -t54 * t314 + (t81 * t343 + t48 * t370) * t173;
t26 = -t53 * t314 + (t81 * t337 + t48 * t371) * t173;
t25 = -t52 * t315 + (t80 * t345 + t47 * t372) * t172;
t24 = -t51 * t315 + (t80 * t339 + t47 * t373) * t172;
t23 = -t50 * t316 + (t79 * t347 + t46 * t374) * t171;
t22 = -t49 * t316 + (t79 * t341 + t46 * t375) * t171;
t18 = t57 - t45 / 0.2e1 - t44 / 0.2e1 - t75 / 0.2e1;
t17 = t56 - t43 / 0.2e1 - t42 / 0.2e1 - t74 / 0.2e1;
t16 = t55 - t41 / 0.2e1 - t40 / 0.2e1 - t73 / 0.2e1;
t15 = (-pkin(3) * t381 + (-t57 + t439 + (-t186 - t406) * t60) * t60) * t190;
t14 = (-pkin(3) * t384 + (-t56 + t440 + (-t186 - t407) * t59) * t59) * t189;
t13 = (-pkin(3) * t387 + (-t55 + t441 + (-t186 - t408) * t58) * t58) * t188;
t12 = ((t18 * t431 + ((qJ(2,1) * t422) - t270 * t152 - (2 * t265) + t276) * t432 + (t119 + t430) * t30) * t60 + ((-0.2e1 * t18 * t185 + t207 * t376 + (-pkin(2) * t182 - t400) * t60) * t60 - (-t207 * t388 / 0.2e1 + t87 * t119) * t173 * t87) * pkin(3)) * t190;
t11 = ((t17 * t431 + ((qJ(2,2) * t422) - t270 * t151 - (2 * t264) + t276) * t433 + (t118 + t430) * t29) * t59 + ((-0.2e1 * t17 * t184 + t206 * t377 + (-pkin(2) * t181 - t401) * t59) * t59 - (-t206 * t389 / 0.2e1 + t86 * t118) * t172 * t86) * pkin(3)) * t189;
t10 = ((t16 * t431 + ((qJ(2,3) * t422) - t270 * t150 - (2 * t263) + t276) * t434 + (t117 + t430) * t28) * t58 + ((-0.2e1 * t16 * t183 + t205 * t378 + (-pkin(2) * t180 - t402) * t58) * t58 - (-t205 * t390 / 0.2e1 + t85 * t117) * t171 * t85) * pkin(3)) * t188;
t9 = -t81 * t15 + t146 * t311 + ((t57 - 0.2e1 * t45 - 0.2e1 * t44 - 0.2e1 * t75) * t429 + (t138 * t149 + 0.2e1 * t349 + (t331 * rSges(3,1) + t328 * rSges(3,2)) * t419) * t60) * t432 + m(3) * ((-rSges(3,1) * t122 + t286 * rSges(3,2)) * t185 + t179 * (t286 * rSges(3,1) + rSges(3,2) * t122));
t8 = -t80 * t14 + t146 * t312 + ((t56 - 0.2e1 * t43 - 0.2e1 * t42 - 0.2e1 * t74) * t428 + (t138 * t148 + 0.2e1 * t350 + (t332 * rSges(3,1) + t329 * rSges(3,2)) * t419) * t59) * t433 + m(3) * ((-rSges(3,1) * t121 + t287 * rSges(3,2)) * t184 + t178 * (t287 * rSges(3,1) + rSges(3,2) * t121));
t7 = -t79 * t13 + t146 * t313 + ((t55 - 0.2e1 * t41 - 0.2e1 * t40 - 0.2e1 * t73) * t427 + (t138 * t147 + 0.2e1 * t351 + (t333 * rSges(3,1) + t330 * rSges(3,2)) * t419) * t58) * t434 + m(3) * ((-rSges(3,1) * t120 + t288 * rSges(3,2)) * t183 + t177 * (t288 * rSges(3,1) + rSges(3,2) * t120));
t6 = t72 * t15 - (t132 * t60 - t87 * t429) * t60 + (-g(3) * t241 - t125 * t235 - t12) * t253;
t5 = t71 * t14 - (t131 * t59 - t86 * t428) * t59 + (-g(3) * t239 - t124 * t233 - t11) * t253;
t4 = t70 * t13 - (t130 * t58 - t85 * t427) * t58 + (-g(3) * t237 - t123 * t231 - t10) * t253;
t3 = -t48 * t15 + t72 * t12 + t81 * t311 + ((-t125 * t203 + t255) * m(3) + t145 * t125 + t156 + t283 * g(3)) * t241 + t235 * ((-t145 + t409) * g(3) + (t283 + t302) * t125) + (-t293 * t185 * t87 - t138 * t388 - t142 * t376) * t87 + (t132 * t439 - t307 * t376 + (t306 * t185 - 0.2e1 * t349 + (-rSges(3,1) * t176 - rSges(3,2) * t182 - t438) * t419) * t87) * t60;
t2 = -t47 * t14 + t71 * t11 + t80 * t312 + ((-t124 * t202 + t255) * m(3) + t144 * t124 + t156 + t284 * g(3)) * t239 + t233 * ((-t144 + t410) * g(3) + (t284 + t302) * t124) + (-t294 * t184 * t86 - t138 * t389 - t141 * t377) * t86 + (t131 * t440 - t307 * t377 + (t306 * t184 - 0.2e1 * t350 + (-rSges(3,1) * t175 - rSges(3,2) * t181 - t437) * t419) * t86) * t59;
t1 = -t46 * t13 + t70 * t10 + t79 * t313 + ((-t123 * t201 + t255) * m(3) + t143 * t123 + t156 + t285 * g(3)) * t237 + t231 * ((-t143 + t411) * g(3) + (t285 + t302) * t123) + (-t295 * t183 * t85 - t138 * t390 - t140 * t378) * t85 + (t130 * t441 - t307 * t378 + (t306 * t183 - 0.2e1 * t351 + (-rSges(3,1) * t174 - rSges(3,2) * t180 - t436) * t419) * t85) * t58;
t19 = [(-g(1) + t229) * m(4) + ((t27 * t379 + t36 * t391) * t229 + (t27 * t380 + t36 * t392) * t228 + (-t235 * t27 + t36 * t90) * t227 + t3 * t379 + t6 * t391) * t190 + ((t196 * t228 * t66 + (t229 * t66 + t9) * t193) * t173 + (t195 * t228 * t64 + (t229 * t64 + t8) * t192) * t172 + (t194 * t228 * t62 + (t229 * t62 + t7) * t191) * t171) * t271 + ((t25 * t382 + t34 * t393) * t229 + (t25 * t383 + t34 * t394) * t228 + (-t233 * t25 + t34 * t89) * t227 + t2 * t382 + t5 * t393) * t189 + ((t23 * t385 + t32 * t395) * t229 + (t23 * t386 + t32 * t396) * t228 + (-t23 * t231 + t32 * t88) * t227 + t1 * t385 + t4 * t395) * t188; (-g(2) + t228) * m(4) + ((t26 * t379 + t35 * t391) * t229 + (t26 * t380 + t35 * t392) * t228 + (-t235 * t26 + t35 * t90) * t227 + t3 * t380 + t6 * t392) * t190 + ((t193 * t229 * t65 + (t228 * t65 + t9) * t196) * t173 + (t192 * t229 * t63 + (t228 * t63 + t8) * t195) * t172 + (t191 * t229 * t61 + (t228 * t61 + t7) * t194) * t171) * t271 + ((t24 * t382 + t33 * t393) * t229 + (t24 * t383 + t33 * t394) * t228 + (-t233 * t24 + t33 * t89) * t227 + t2 * t383 + t5 * t394) * t189 + ((t22 * t385 + t31 * t395) * t229 + (t22 * t386 + t31 * t396) * t228 + (-t22 * t231 + t31 * t88) * t227 + t1 * t386 + t4 * t396) * t188; (-g(3) + t227) * m(4) + ((-t235 * t39 + t69 * t90) * t227 - t235 * t3 + t90 * t6 + (t228 * t53 + t229 * t54) * t69 * t112 + ((-t193 * t308 + t39 * t99) * t229 + (-t196 * t308 + t39 * t98) * t228) * t173) * t190 + ((-t233 * t38 + t68 * t89) * t227 - t233 * t2 + t89 * t5 + (t228 * t51 + t229 * t52) * t68 * t111 + ((-t192 * t309 + t38 * t97) * t229 + (-t195 * t309 + t38 * t96) * t228) * t172) * t189 + ((-t231 * t37 + t67 * t88) * t227 - t231 * t1 + t88 * t4 + (t228 * t49 + t229 * t50) * t67 * t110 + ((-t191 * t310 + t37 * t95) * t229 + (-t194 * t310 + t37 * t94) * t228) * t171) * t188;];
tauX  = t19;
