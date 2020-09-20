% Calculate vector of inverse dynamics forces for parallel robot
% P3RRRRR2G2P2A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2020-03-09 21:10
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRRRR2G2P2A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G2P2A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRRRR2G2P2A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRRRR2G2P2A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G2P2A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR2G2P2A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G2P2A0_invdyn_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR2G2P2A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRRRR2G2P2A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRRRR2G2P2A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G2P2A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G2P2A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:08:30
% EndTime: 2020-03-09 21:08:37
% DurationCPUTime: 7.22s
% Computational Cost: add. (16026->502), mult. (32016->913), div. (9969->14), fcn. (29079->60), ass. (0->388)
t267 = xDP(2);
t255 = cos(qJ(3,1));
t224 = 0.1e1 / t255 ^ 2;
t247 = sin(qJ(2,1));
t215 = 0.1e1 / t247;
t273 = 0.1e1 / pkin(2);
t275 = 0.1e1 / pkin(1);
t343 = t273 * t275;
t299 = t215 * t343;
t284 = t224 * t299;
t236 = legFrame(1,2);
t201 = sin(t236);
t204 = cos(t236);
t248 = sin(qJ(1,1));
t256 = cos(qJ(2,1));
t257 = cos(qJ(1,1));
t138 = t247 * t257 + t248 * t256;
t222 = t255 ^ 2;
t333 = pkin(2) * t138 * t222;
t246 = sin(qJ(3,1));
t439 = pkin(1) * t256;
t336 = t246 * t439;
t361 = t204 * t246;
t442 = pkin(1) * t248;
t80 = t201 * t333 + (-pkin(2) * t361 + t201 * t442) * t255 - t204 * t336;
t74 = t80 * t267 * t284;
t268 = xDP(1);
t364 = t201 * t246;
t81 = -t204 * t333 + (-pkin(2) * t364 - t204 * t442) * t255 - t201 * t336;
t75 = t81 * t268 * t284;
t266 = xDP(3);
t298 = t266 * t343;
t231 = qJ(2,1) + qJ(3,1);
t197 = qJ(1,1) + t231;
t169 = cos(t197);
t232 = qJ(2,1) - qJ(3,1);
t198 = qJ(1,1) + t232;
t170 = cos(t198);
t181 = sin(t231);
t182 = sin(t232);
t458 = -0.2e1 * pkin(1);
t394 = (t257 * t458 + (-t169 - t170) * pkin(2)) / (t181 + t182);
t93 = t298 * t394;
t48 = t75 + t74 + t93;
t388 = t138 * t255;
t101 = -t201 * t388 + t361;
t102 = t204 * t388 + t364;
t223 = 0.1e1 / t255;
t358 = t215 * t275;
t300 = t223 * t358;
t344 = t266 * t275;
t233 = qJ(1,1) + qJ(2,1);
t192 = cos(t233);
t368 = t192 * t215;
t66 = t344 * t368 + (t101 * t267 + t102 * t268) * t300;
t433 = -t48 - t66;
t471 = t433 ^ 2;
t252 = cos(qJ(3,2));
t221 = 0.1e1 / t252 ^ 2;
t244 = sin(qJ(2,2));
t214 = 0.1e1 / t244;
t301 = t214 * t343;
t285 = t221 * t301;
t235 = legFrame(2,2);
t200 = sin(t235);
t203 = cos(t235);
t245 = sin(qJ(1,2));
t253 = cos(qJ(2,2));
t254 = cos(qJ(1,2));
t137 = t244 * t254 + t245 * t253;
t219 = t252 ^ 2;
t334 = pkin(2) * t137 * t219;
t243 = sin(qJ(3,2));
t440 = pkin(1) * t253;
t337 = t243 * t440;
t362 = t203 * t243;
t443 = pkin(1) * t245;
t78 = t200 * t334 + (-pkin(2) * t362 + t200 * t443) * t252 - t203 * t337;
t72 = t78 * t267 * t285;
t365 = t200 * t243;
t79 = -t203 * t334 + (-pkin(2) * t365 - t203 * t443) * t252 - t200 * t337;
t73 = t79 * t268 * t285;
t228 = qJ(2,2) + qJ(3,2);
t195 = qJ(1,2) + t228;
t167 = cos(t195);
t229 = qJ(2,2) - qJ(3,2);
t196 = qJ(1,2) + t229;
t168 = cos(t196);
t178 = sin(t228);
t179 = sin(t229);
t395 = (t254 * t458 + (-t167 - t168) * pkin(2)) / (t178 + t179);
t92 = t298 * t395;
t47 = t73 + t72 + t92;
t389 = t137 * t252;
t100 = t203 * t389 + t365;
t220 = 0.1e1 / t252;
t359 = t214 * t275;
t302 = t220 * t359;
t230 = qJ(1,2) + qJ(2,2);
t189 = cos(t230);
t370 = t189 * t214;
t99 = -t200 * t389 + t362;
t65 = t344 * t370 + (t100 * t268 + t267 * t99) * t302;
t434 = -t47 - t65;
t470 = t434 ^ 2;
t249 = cos(qJ(3,3));
t218 = 0.1e1 / t249 ^ 2;
t241 = sin(qJ(2,3));
t213 = 0.1e1 / t241;
t303 = t213 * t343;
t286 = t218 * t303;
t234 = legFrame(3,2);
t199 = sin(t234);
t202 = cos(t234);
t242 = sin(qJ(1,3));
t250 = cos(qJ(2,3));
t251 = cos(qJ(1,3));
t136 = t241 * t251 + t242 * t250;
t216 = t249 ^ 2;
t335 = pkin(2) * t136 * t216;
t240 = sin(qJ(3,3));
t441 = pkin(1) * t250;
t338 = t240 * t441;
t363 = t202 * t240;
t444 = pkin(1) * t242;
t76 = t199 * t335 + (-pkin(2) * t363 + t199 * t444) * t249 - t202 * t338;
t70 = t76 * t267 * t286;
t366 = t199 * t240;
t77 = -t202 * t335 + (-pkin(2) * t366 - t202 * t444) * t249 - t199 * t338;
t71 = t77 * t268 * t286;
t225 = qJ(2,3) + qJ(3,3);
t193 = qJ(1,3) + t225;
t165 = cos(t193);
t226 = qJ(2,3) - qJ(3,3);
t194 = qJ(1,3) + t226;
t166 = cos(t194);
t175 = sin(t225);
t176 = sin(t226);
t396 = (t251 * t458 + (-t165 - t166) * pkin(2)) / (t175 + t176);
t91 = t298 * t396;
t46 = t71 + t70 + t91;
t217 = 0.1e1 / t249;
t360 = t213 * t275;
t304 = t217 * t360;
t227 = qJ(1,3) + qJ(2,3);
t186 = cos(t227);
t372 = t186 * t213;
t390 = t136 * t249;
t97 = -t199 * t390 + t363;
t98 = t202 * t390 + t366;
t64 = t344 * t372 + (t267 * t97 + t268 * t98) * t304;
t435 = -t46 - t64;
t469 = t435 ^ 2;
t454 = m(3) * pkin(1);
t468 = m(1) * rSges(1,1) + m(2) * pkin(1) + t454;
t341 = 0.2e1 * pkin(1);
t451 = m(3) * rSges(3,2);
t174 = rSges(3,1) * t451 - Icges(3,4);
t271 = 0.2e1 * qJ(3,1);
t208 = sin(t271);
t211 = cos(t271);
t456 = rSges(3,2) ^ 2;
t457 = rSges(3,1) ^ 2;
t157 = (-t456 + t457) * m(3) - Icges(3,1) + Icges(3,2);
t445 = t157 / 0.2e1;
t467 = -t174 * t208 + t211 * t445;
t270 = 0.2e1 * qJ(3,2);
t207 = sin(t270);
t210 = cos(t270);
t466 = -t174 * t207 + t210 * t445;
t269 = 0.2e1 * qJ(3,3);
t206 = sin(t269);
t209 = cos(t269);
t465 = -t174 * t206 + t209 * t445;
t450 = m(3) * rSges(3,3);
t164 = m(2) * rSges(2,2) - t450;
t265 = m(2) * rSges(2,1);
t461 = -rSges(3,1) * t255 + rSges(3,2) * t246;
t464 = -t164 * t247 - (m(3) * t461 - t265) * t256;
t460 = -rSges(3,1) * t252 + rSges(3,2) * t243;
t463 = -t164 * t244 - (m(3) * t460 - t265) * t253;
t459 = -rSges(3,1) * t249 + rSges(3,2) * t240;
t462 = -t164 * t241 - (m(3) * t459 - t265) * t250;
t61 = t64 ^ 2;
t62 = t65 ^ 2;
t63 = t66 ^ 2;
t455 = -0.2e1 * t157;
t453 = m(1) * rSges(1,2);
t452 = m(3) * rSges(3,1);
t449 = rSges(3,3) * g(3);
t448 = pkin(1) * t61;
t447 = pkin(1) * t62;
t446 = pkin(1) * t63;
t25 = t71 / 0.2e1 + t70 / 0.2e1 + t91 / 0.2e1 + t64;
t438 = t46 * t25;
t26 = t73 / 0.2e1 + t72 / 0.2e1 + t92 / 0.2e1 + t65;
t437 = t47 * t26;
t27 = t75 / 0.2e1 + t74 / 0.2e1 + t93 / 0.2e1 + t66;
t436 = t48 * t27;
t139 = g(1) * t202 - g(2) * t199;
t432 = rSges(3,1) * t139;
t140 = g(1) * t203 - g(2) * t200;
t431 = rSges(3,1) * t140;
t141 = g(1) * t204 - g(2) * t201;
t430 = rSges(3,1) * t141;
t423 = rSges(3,3) * t139;
t422 = rSges(3,3) * t140;
t421 = rSges(3,3) * t141;
t357 = t217 * t273;
t106 = (t199 * t268 + t202 * t267) * t357;
t420 = t106 * t435;
t355 = t220 * t273;
t107 = (t200 * t268 + t203 * t267) * t355;
t419 = t107 * t434;
t353 = t223 * t273;
t108 = (t201 * t268 + t204 * t267) * t353;
t418 = t108 * t433;
t103 = t106 ^ 2;
t417 = (t103 / 0.2e1 + t438) * t241;
t104 = t107 ^ 2;
t416 = (t104 / 0.2e1 + t437) * t244;
t105 = t108 ^ 2;
t415 = (t105 / 0.2e1 + t436) * t247;
t414 = t216 * t435;
t172 = -rSges(3,2) * t450 + Icges(3,6);
t173 = rSges(3,1) * t450 - Icges(3,5);
t118 = t172 * t249 - t173 * t240;
t316 = t273 * t396;
t293 = pkin(1) * t241 + rSges(3,3);
t94 = (-t293 * t451 + Icges(3,6)) * t249 - t240 * (t293 * t452 - Icges(3,5));
t413 = t217 * (t118 * t316 + t94 * t372) * t275;
t412 = t217 * t97;
t411 = t217 * t98;
t410 = t219 * t434;
t119 = t172 * t252 - t173 * t243;
t315 = t273 * t395;
t292 = pkin(1) * t244 + rSges(3,3);
t95 = (-t292 * t451 + Icges(3,6)) * t252 - t243 * (t292 * t452 - Icges(3,5));
t409 = t220 * (t119 * t315 + t95 * t370) * t275;
t408 = t220 * t99;
t407 = t222 * t433;
t120 = t172 * t255 - t173 * t246;
t314 = t273 * t394;
t291 = pkin(1) * t247 + rSges(3,3);
t96 = (-t291 * t451 + Icges(3,6)) * t255 - t246 * (t291 * t452 - Icges(3,5));
t406 = t223 * (t120 * t314 + t96 * t368) * t275;
t405 = t100 * t220;
t404 = t101 * t223;
t403 = t102 * t223;
t402 = t103 * t217;
t401 = t104 * t220;
t400 = t105 * t223;
t399 = t106 * t250;
t398 = t107 * t253;
t397 = t108 * t256;
t387 = t157 * t206;
t386 = t157 * t207;
t385 = t157 * t208;
t381 = t173 * t249;
t380 = t173 * t252;
t379 = t173 * t255;
t375 = t174 * t209;
t374 = t174 * t210;
t373 = t174 * t211;
t237 = xDDP(3);
t371 = t186 * t237;
t369 = t189 * t237;
t367 = t192 * t237;
t356 = t218 * t273;
t354 = t221 * t273;
t352 = t224 * t273;
t351 = t237 * t273;
t350 = t240 * t241;
t349 = t243 * t244;
t348 = t246 * t247;
t347 = t249 * t250;
t346 = t252 * t253;
t345 = t255 * t256;
t342 = t468 * g(3);
t340 = m(3) * t449;
t339 = t454 / 0.2e1;
t332 = t435 * t399;
t331 = t434 * t398;
t330 = t433 * t397;
t329 = t76 * t356;
t328 = t77 * t356;
t320 = t456 + t457;
t277 = Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1 + (0.2e1 * rSges(3,3) ^ 2 + t320) * m(3) / 0.2e1;
t82 = t277 + t465;
t327 = t82 * t356;
t326 = t78 * t354;
t325 = t79 * t354;
t83 = t277 + t466;
t324 = t83 * t354;
t323 = t80 * t352;
t322 = t81 * t352;
t84 = t277 + t467;
t321 = t84 * t352;
t319 = t106 * t350;
t318 = t107 * t349;
t317 = t108 * t348;
t313 = t118 * t356;
t312 = t119 * t354;
t311 = t120 * t352;
t310 = t199 * t357;
t309 = t200 * t355;
t308 = t201 * t353;
t307 = t202 * t357;
t306 = t203 * t355;
t305 = t204 * t353;
t290 = -pkin(1) * t452 / 0.2e1;
t280 = t164 * t250 + t241 * t265;
t279 = t164 * t253 + t244 * t265;
t278 = t164 * t256 + t247 * t265;
t276 = Icges(1,3) + (m(3) + m(2)) * pkin(1) ^ 2 + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + t277;
t272 = pkin(2) ^ 2;
t263 = rSges(2,1) * g(3);
t262 = rSges(3,1) * g(3);
t261 = rSges(2,2) * g(3);
t260 = rSges(3,2) * g(3);
t239 = xDDP(1);
t238 = xDDP(2);
t205 = g(3) * t453;
t191 = cos(t232);
t190 = cos(t231);
t188 = cos(t229);
t187 = cos(t228);
t185 = cos(t226);
t184 = cos(t225);
t183 = sin(t233);
t180 = sin(t230);
t177 = sin(t227);
t162 = t320 * m(3) + Icges(3,3);
t129 = rSges(3,2) * t141;
t128 = rSges(3,2) * t140;
t127 = rSges(3,2) * t139;
t117 = m(2) * (rSges(2,2) * t141 + t263);
t116 = m(2) * (rSges(2,2) * t140 + t263);
t115 = m(2) * (rSges(2,2) * t139 + t263);
t114 = m(2) * (-rSges(2,1) * t141 + t261);
t113 = m(2) * (-rSges(2,1) * t140 + t261);
t112 = m(2) * (-rSges(2,1) * t139 + t261);
t69 = t464 * pkin(1) + t84;
t68 = t463 * pkin(1) + t83;
t67 = t462 * pkin(1) + t82;
t60 = t464 * t341 + t276 + t467;
t59 = t463 * t341 + t276 + t466;
t58 = t462 * t341 + t276 + t465;
t54 = t162 * t308 + (t81 * t311 + t96 * t403) * t358;
t53 = t162 * t305 + (t80 * t311 + t96 * t404) * t358;
t52 = t162 * t309 + (t79 * t312 + t95 * t405) * t359;
t51 = t162 * t306 + (t78 * t312 + t95 * t408) * t359;
t50 = t162 * t310 + (t77 * t313 + t94 * t411) * t360;
t49 = t162 * t307 + (t76 * t313 + t94 * t412) * t360;
t45 = (t84 * t314 + t69 * t368) * t275;
t44 = (t83 * t315 + t68 * t370) * t275;
t43 = (t82 * t316 + t67 * t372) * t275;
t42 = (t69 * t314 + t60 * t368) * t275;
t41 = (t68 * t315 + t59 * t370) * t275;
t40 = (t67 * t316 + t58 * t372) * t275;
t39 = t120 * t308 + (t81 * t321 + t69 * t403) * t358;
t38 = t120 * t305 + (t80 * t321 + t69 * t404) * t358;
t37 = t119 * t309 + (t79 * t324 + t68 * t405) * t359;
t36 = t119 * t306 + (t78 * t324 + t68 * t408) * t359;
t35 = t118 * t310 + (t77 * t327 + t67 * t411) * t360;
t34 = t118 * t307 + (t76 * t327 + t67 * t412) * t360;
t24 = t96 * t308 + (t69 * t322 + t60 * t403) * t358;
t23 = t96 * t305 + (t69 * t323 + t60 * t404) * t358;
t22 = t95 * t309 + (t68 * t325 + t59 * t405) * t359;
t21 = t95 * t306 + (t68 * t326 + t59 * t408) * t359;
t20 = t94 * t310 + (t67 * t328 + t58 * t411) * t360;
t19 = t94 * t307 + (t67 * t329 + t58 * t412) * t360;
t15 = (-t63 * t439 + (-t255 * t471 - t400) * pkin(2)) * t358;
t14 = (-t62 * t440 + (-t252 * t470 - t401) * pkin(2)) * t359;
t13 = (-t61 * t441 + (-t249 * t469 - t402) * pkin(2)) * t360;
t12 = (-t272 * t407 + (pkin(1) * t66 + (0.2e1 * t27 * t345 - t317) * pkin(2)) * pkin(1)) * t223 * t66 * t299 + (-pkin(2) * t407 + (-t345 * t433 - t317) * pkin(1)) * t48 * t300 + ((pkin(1) * t348 * t433 + pkin(2) * t108) * t255 + pkin(1) * t397) * t224 * t108 * t358;
t11 = (-t272 * t410 + (pkin(1) * t65 + (0.2e1 * t26 * t346 - t318) * pkin(2)) * pkin(1)) * t220 * t65 * t301 + (-pkin(2) * t410 + (-t346 * t434 - t318) * pkin(1)) * t47 * t302 + ((pkin(1) * t349 * t434 + pkin(2) * t107) * t252 + pkin(1) * t398) * t221 * t107 * t359;
t10 = (-t272 * t414 + (pkin(1) * t64 + (0.2e1 * t25 * t347 - t319) * pkin(2)) * pkin(1)) * t217 * t64 * t303 + (-pkin(2) * t414 + (-t347 * t435 - t319) * pkin(1)) * t46 * t304 + ((pkin(1) * t350 * t435 + pkin(2) * t106) * t249 + pkin(1) * t399) * t218 * t106 * t360;
t9 = -t96 * t15 - t120 * t12 + t162 * t246 * t400 + ((g(1) * t201 + g(2) * t204) * t461 + (g(3) * t192 + t141 * t183) * (rSges(3,1) * t246 + rSges(3,2) * t255)) * m(3) + (t385 / 0.2e1 + t373) * t471 + (t182 * t290 + (rSges(3,1) * t181 + (t190 + t191) * rSges(3,2)) * t339) * t63;
t8 = -t95 * t14 - t119 * t11 + t162 * t243 * t401 + ((g(1) * t200 + g(2) * t203) * t460 + (g(3) * t189 + t140 * t180) * (rSges(3,1) * t243 + rSges(3,2) * t252)) * m(3) + (t386 / 0.2e1 + t374) * t470 + (t179 * t290 + (rSges(3,1) * t178 + (t187 + t188) * rSges(3,2)) * t339) * t62;
t7 = -t94 * t13 - t118 * t10 + t162 * t240 * t402 + ((g(1) * t199 + g(2) * t202) * t459 + (g(3) * t186 + t139 * t177) * (rSges(3,1) * t240 + rSges(3,2) * t249)) * m(3) + (t387 / 0.2e1 + t375) * t469 + (t176 * t290 + (rSges(3,1) * t175 + (t184 + t185) * rSges(3,2)) * t339) * t61;
t6 = t114 * t192 + t183 * t117 - t84 * t12 - t69 * t15 - (-0.2e1 * t373 - t385) * t418 + t278 * t446 + (-t379 + (t120 * t223 - t172) * t246) * t105 + ((t461 * t141 - t449) * t192 + t183 * (-t461 * g(3) - t421) + ((-t191 / 0.2e1 + t190 / 0.2e1) * rSges(3,2) + (t182 / 0.2e1 + t181 / 0.2e1) * rSges(3,1)) * t446) * m(3);
t5 = -t83 * t11 + t113 * t189 + t180 * t116 - t68 * t14 - (-0.2e1 * t374 - t386) * t419 + t279 * t447 + (-t380 + (t119 * t220 - t172) * t243) * t104 + ((t460 * t140 - t449) * t189 + t180 * (-t460 * g(3) - t422) + ((-t188 / 0.2e1 + t187 / 0.2e1) * rSges(3,2) + (t179 / 0.2e1 + t178 / 0.2e1) * rSges(3,1)) * t447) * m(3);
t4 = -t82 * t10 + t112 * t186 + t177 * t115 - t67 * t13 - (-0.2e1 * t375 - t387) * t420 + t280 * t448 + (-t381 + (t118 * t217 - t172) * t240) * t103 + ((t459 * t139 - t449) * t186 + t177 * (-t459 * g(3) - t423) + ((-t185 / 0.2e1 + t184 / 0.2e1) * rSges(3,2) + (t176 / 0.2e1 + t175 / 0.2e1) * rSges(3,1)) * t448) * m(3);
t3 = -t60 * t15 - t69 * t12 + (t141 * t453 + t342) * t248 - (t246 * t255 * t455 + (-0.4e1 * t222 + 0.2e1) * t174) * t418 + (-t379 + (t223 * t96 - t172) * t246) * t105 - t278 * t341 * t436 + (-(t260 + t430) * t170 / 0.2e1 - (t129 - t262) * sin(t198) / 0.2e1 + (t260 - t430) * t169 / 0.2e1 + (t129 + t262) * sin(t197) / 0.2e1 + ((-rSges(3,1) * t415 + rSges(3,2) * t330) * t255 + (rSges(3,1) * t330 + rSges(3,2) * t415) * t246) * t341) * m(3) + (t114 - t340) * t192 + t183 * (-m(3) * t421 + t117) + (-t141 * t468 + t205) * t257;
t2 = -t59 * t14 - t68 * t11 + (t140 * t453 + t342) * t245 - (t243 * t252 * t455 + (-0.4e1 * t219 + 0.2e1) * t174) * t419 + (-t380 + (t220 * t95 - t172) * t243) * t104 - t279 * t341 * t437 + (-(t260 + t431) * t168 / 0.2e1 - (t128 - t262) * sin(t196) / 0.2e1 + (t260 - t431) * t167 / 0.2e1 + (t128 + t262) * sin(t195) / 0.2e1 + ((-rSges(3,1) * t416 + rSges(3,2) * t331) * t252 + (rSges(3,1) * t331 + rSges(3,2) * t416) * t243) * t341) * m(3) + (t113 - t340) * t189 + t180 * (-m(3) * t422 + t116) + (-t140 * t468 + t205) * t254;
t1 = -t58 * t13 - t67 * t10 + (t139 * t453 + t342) * t242 - (t240 * t249 * t455 + (-0.4e1 * t216 + 0.2e1) * t174) * t420 + (-t381 + (t217 * t94 - t172) * t240) * t103 - t280 * t341 * t438 + (-(t260 + t432) * t166 / 0.2e1 - (t127 - t262) * sin(t194) / 0.2e1 + (t260 - t432) * t165 / 0.2e1 + (t127 + t262) * sin(t193) / 0.2e1 + ((-rSges(3,1) * t417 + rSges(3,2) * t332) * t249 + (rSges(3,1) * t332 + rSges(3,2) * t417) * t240) * t341) * m(3) + (t112 - t340) * t186 + t177 * (-m(3) * t423 + t115) + (-t139 * t468 + t205) * t251;
t16 = [(-g(1) + t239) * m(4) + ((t204 * t238 * t54 + (t239 * t54 + t9) * t201) * t223 + (t203 * t238 * t52 + (t239 * t52 + t8) * t200) * t220 + (t202 * t238 * t50 + (t239 * t50 + t7) * t199) * t217) * t273 + ((t35 * t396 + t37 * t395 + t39 * t394) * t351 + ((t24 * t403 + t39 * t322) * t239 + (t24 * t404 + t39 * t323) * t238 + t24 * t367 + t3 * t403 + t6 * t322) * t215 + ((t22 * t405 + t37 * t325) * t239 + (t22 * t408 + t37 * t326) * t238 + t22 * t369 + t2 * t405 + t5 * t325) * t214 + ((t20 * t411 + t35 * t328) * t239 + (t20 * t412 + t35 * t329) * t238 + t20 * t371 + t1 * t411 + t4 * t328) * t213) * t275; (-g(2) + t238) * m(4) + ((t201 * t239 * t53 + (t238 * t53 + t9) * t204) * t223 + (t200 * t239 * t51 + (t238 * t51 + t8) * t203) * t220 + (t199 * t239 * t49 + (t238 * t49 + t7) * t202) * t217) * t273 + ((t34 * t396 + t36 * t395 + t38 * t394) * t351 + ((t23 * t403 + t38 * t322) * t239 + (t23 * t404 + t38 * t323) * t238 + t23 * t367 + t3 * t404 + t6 * t323) * t215 + ((t21 * t405 + t36 * t325) * t239 + (t21 * t408 + t36 * t326) * t238 + t21 * t369 + t2 * t408 + t5 * t326) * t214 + ((t19 * t411 + t34 * t328) * t239 + (t19 * t412 + t34 * t329) * t238 + t19 * t371 + t1 * t412 + t4 * t329) * t213) * t275; (-g(3) + t237) * m(4) + ((t199 * t413 + t200 * t409 + t201 * t406) * t239 + (t202 * t413 + t203 * t409 + t204 * t406) * t238) * t273 + ((t4 * t396 + t5 * t395 + t6 * t394 + (t45 * t394 + t44 * t395 + t43 * t396) * t237) * t273 + ((t45 * t322 + t42 * t403) * t239 + (t45 * t323 + t42 * t404) * t238 + (t237 * t42 + t3) * t192) * t215 + ((t44 * t325 + t41 * t405) * t239 + (t44 * t326 + t41 * t408) * t238 + (t237 * t41 + t2) * t189) * t214 + ((t43 * t328 + t40 * t411) * t239 + (t43 * t329 + t40 * t412) * t238 + (t237 * t40 + t1) * t186) * t213) * t275;];
tauX  = t16;
