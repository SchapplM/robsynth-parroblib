% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRRR8V2G4A0
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
% Datum: 2020-08-06 18:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRRR8V2G4A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G4A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V2G4A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR8V2G4A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G4A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V2G4A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G4A0_invdyn_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G4A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V2G4A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR8V2G4A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G4A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G4A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:14:07
% EndTime: 2020-08-06 18:14:21
% DurationCPUTime: 13.94s
% Computational Cost: add. (63081->599), mult. (129003->1129), div. (5112->10), fcn. (147594->34), ass. (0->442)
t277 = cos(qJ(2,1));
t284 = pkin(7) + pkin(6);
t218 = t277 * t284;
t271 = sin(qJ(2,1));
t187 = pkin(2) * t271 - t218;
t251 = sin(pkin(4));
t253 = cos(pkin(4));
t270 = sin(qJ(3,1));
t360 = t253 * t270;
t153 = pkin(3) * t360 + t187 * t251;
t276 = cos(qJ(3,1));
t371 = t251 * t271;
t249 = t276 ^ 2;
t465 = pkin(3) * t249;
t123 = 0.1e1 / (pkin(2) * t360 + t153 * t276 + t371 * t465);
t275 = cos(qJ(2,2));
t217 = t275 * t284;
t269 = sin(qJ(2,2));
t186 = pkin(2) * t269 - t217;
t268 = sin(qJ(3,2));
t362 = t253 * t268;
t152 = pkin(3) * t362 + t186 * t251;
t274 = cos(qJ(3,2));
t373 = t251 * t269;
t248 = t274 ^ 2;
t466 = pkin(3) * t248;
t122 = 0.1e1 / (pkin(2) * t362 + t152 * t274 + t373 * t466);
t273 = cos(qJ(2,3));
t216 = t273 * t284;
t267 = sin(qJ(2,3));
t185 = pkin(2) * t267 - t216;
t266 = sin(qJ(3,3));
t364 = t253 * t266;
t151 = pkin(3) * t364 + t185 * t251;
t272 = cos(qJ(3,3));
t375 = t251 * t267;
t247 = t272 ^ 2;
t467 = pkin(3) * t247;
t121 = 0.1e1 / (pkin(2) * t364 + t151 * t272 + t375 * t467);
t496 = -rSges(3,1) * t276 + rSges(3,2) * t270;
t495 = -rSges(3,1) * t274 + rSges(3,2) * t268;
t494 = -rSges(3,1) * t272 + rSges(3,2) * t266;
t259 = legFrame(1,1);
t228 = sin(t259);
t234 = cos(t259);
t262 = legFrame(1,2);
t237 = sin(t262);
t240 = cos(t262);
t150 = g(1) * t237 + (-g(2) * t228 + g(3) * t234) * t240;
t256 = legFrame(1,3);
t225 = sin(t256);
t231 = cos(t256);
t383 = t234 * t237;
t386 = t228 * t237;
t462 = g(1) * t240;
t110 = -t225 * t462 + (-t225 * t386 + t231 * t234) * g(2) + (t225 * t383 + t228 * t231) * g(3);
t111 = t231 * t462 + (t225 * t234 + t231 * t386) * g(2) + (t225 * t228 - t231 * t383) * g(3);
t250 = sin(pkin(8));
t252 = cos(pkin(8));
t301 = t110 * t252 - t111 * t250;
t493 = -t150 * t253 + t301 * t251;
t258 = legFrame(2,1);
t227 = sin(t258);
t233 = cos(t258);
t261 = legFrame(2,2);
t236 = sin(t261);
t239 = cos(t261);
t149 = g(1) * t236 + (-g(2) * t227 + g(3) * t233) * t239;
t255 = legFrame(2,3);
t224 = sin(t255);
t230 = cos(t255);
t384 = t233 * t236;
t387 = t227 * t236;
t463 = g(1) * t239;
t108 = -t224 * t463 + (-t224 * t387 + t230 * t233) * g(2) + (t224 * t384 + t227 * t230) * g(3);
t109 = t230 * t463 + (t224 * t233 + t230 * t387) * g(2) + (t224 * t227 - t230 * t384) * g(3);
t303 = t108 * t252 - t109 * t250;
t492 = -t149 * t253 + t303 * t251;
t257 = legFrame(3,1);
t226 = sin(t257);
t232 = cos(t257);
t260 = legFrame(3,2);
t235 = sin(t260);
t238 = cos(t260);
t148 = g(1) * t235 + (-g(2) * t226 + g(3) * t232) * t238;
t254 = legFrame(3,3);
t223 = sin(t254);
t229 = cos(t254);
t385 = t232 * t235;
t388 = t226 * t235;
t464 = g(1) * t238;
t106 = -t223 * t464 + (-t223 * t388 + t229 * t232) * g(2) + (t223 * t385 + t226 * t229) * g(3);
t107 = t229 * t464 + (t223 * t232 + t229 * t388) * g(2) + (t223 * t226 - t229 * t385) * g(3);
t305 = t106 * t252 - t107 * t250;
t491 = -t148 * t253 + t305 * t251;
t490 = t150 * t251 + t301 * t253;
t489 = t149 * t251 + t303 * t253;
t488 = t148 * t251 + t305 * t253;
t482 = m(3) * rSges(3,1);
t345 = rSges(3,2) * t482;
t219 = -Icges(3,4) + t345;
t285 = pkin(2) * m(3);
t341 = t285 / 0.2e1;
t319 = rSges(3,1) * t341;
t487 = t219 * t247 + t266 * t319;
t486 = t219 * t248 + t268 * t319;
t485 = t219 * t249 + t270 * t319;
t484 = 0.2e1 * pkin(2);
t357 = t267 * t272;
t157 = pkin(3) * t357 + t185;
t213 = pkin(3) * t272 + pkin(2);
t181 = t213 * t364;
t370 = t251 * t272;
t136 = 0.1e1 / (t157 * t370 + t181);
t281 = xDP(3);
t282 = xDP(2);
t283 = xDP(1);
t381 = t238 * t283;
t160 = -t223 * t250 + t229 * t252;
t163 = t223 * t252 + t229 * t250;
t115 = t160 * t385 - t163 * t226;
t363 = t253 * t267;
t169 = t250 * t363 - t273 * t252;
t172 = t273 * t250 + t252 * t363;
t124 = -t169 * t223 + t172 * t229;
t299 = t169 * t229 + t172 * t223;
t76 = (t124 * t385 - t226 * t299) * t266 + t115 * t370;
t118 = t160 * t388 + t163 * t232;
t79 = (-t124 * t388 - t299 * t232) * t266 - t118 * t370;
t97 = t160 * t370 + t266 * (t160 * t363 + t163 * t273);
t58 = -t121 * t97 * t381 + (t281 * t76 + t282 * t79) * t136;
t55 = t58 ^ 2;
t355 = t269 * t274;
t158 = pkin(3) * t355 + t186;
t214 = pkin(3) * t274 + pkin(2);
t182 = t214 * t362;
t368 = t251 * t274;
t137 = 0.1e1 / (t158 * t368 + t182);
t379 = t239 * t283;
t161 = -t224 * t250 + t230 * t252;
t164 = t224 * t252 + t230 * t250;
t116 = t161 * t384 - t164 * t227;
t361 = t253 * t269;
t170 = t250 * t361 - t275 * t252;
t173 = t275 * t250 + t252 * t361;
t125 = -t170 * t224 + t173 * t230;
t298 = t170 * t230 + t173 * t224;
t77 = (t125 * t384 - t227 * t298) * t268 + t116 * t368;
t119 = t161 * t387 + t164 * t233;
t80 = (-t125 * t387 - t298 * t233) * t268 - t119 * t368;
t98 = t161 * t368 + t268 * (t161 * t361 + t164 * t275);
t59 = -t122 * t98 * t379 + (t281 * t77 + t282 * t80) * t137;
t56 = t59 ^ 2;
t353 = t271 * t276;
t159 = pkin(3) * t353 + t187;
t215 = pkin(3) * t276 + pkin(2);
t183 = t215 * t360;
t366 = t251 * t276;
t138 = 0.1e1 / (t159 * t366 + t183);
t377 = t240 * t283;
t162 = -t225 * t250 + t231 * t252;
t165 = t225 * t252 + t231 * t250;
t117 = t162 * t383 - t165 * t228;
t359 = t253 * t271;
t171 = t250 * t359 - t277 * t252;
t174 = t277 * t250 + t252 * t359;
t126 = -t171 * t225 + t174 * t231;
t297 = t171 * t231 + t174 * t225;
t78 = (t126 * t383 - t228 * t297) * t270 + t117 * t366;
t120 = t162 * t386 + t165 * t234;
t81 = (-t126 * t386 - t297 * t234) * t270 - t120 * t366;
t99 = t162 * t366 + t270 * (t162 * t359 + t165 * t277);
t60 = -t123 * t99 * t377 + (t281 * t78 + t282 * t81) * t138;
t57 = t60 ^ 2;
t483 = -0.2e1 * t219;
t178 = t213 * t267 - t216;
t356 = t267 * t284;
t395 = (t213 * t273 + t356) * t253;
t100 = -t160 * t395 + t163 * t178;
t139 = 0.1e1 / (t178 * t370 + t181);
t289 = 0.1e1 / pkin(3);
t410 = t139 * t289;
t82 = -t118 * t395 + (-t160 * t232 + t163 * t388) * t178;
t85 = t115 * t395 - (t160 * t226 + t163 * t385) * t178;
t64 = (t100 * t381 + t281 * t85 + t282 * t82) * t410;
t481 = pkin(3) * t64;
t179 = t214 * t269 - t217;
t354 = t269 * t284;
t394 = (t214 * t275 + t354) * t253;
t101 = -t161 * t394 + t164 * t179;
t140 = 0.1e1 / (t179 * t368 + t182);
t409 = t140 * t289;
t83 = -t119 * t394 + (-t161 * t233 + t164 * t387) * t179;
t86 = t116 * t394 - (t161 * t227 + t164 * t384) * t179;
t65 = (t101 * t379 + t281 * t86 + t282 * t83) * t409;
t480 = pkin(3) * t65;
t180 = t215 * t271 - t218;
t352 = t271 * t284;
t393 = (t215 * t277 + t352) * t253;
t102 = -t162 * t393 + t165 * t180;
t141 = 0.1e1 / (t180 * t366 + t183);
t408 = t141 * t289;
t84 = -t120 * t393 + (-t162 * t234 + t165 * t386) * t180;
t87 = t117 * t393 - (t162 * t228 + t165 * t383) * t180;
t66 = (t102 * t377 + t281 * t87 + t282 * t84) * t408;
t479 = pkin(3) * t66;
t286 = rSges(3,2) ^ 2;
t287 = rSges(3,1) ^ 2;
t200 = (-t286 + t287) * m(3) + Icges(3,2) - Icges(3,1);
t478 = t200 / 0.2e1;
t278 = pkin(6) + rSges(3,3);
t471 = m(3) * t278;
t204 = rSges(3,2) * t471 - Icges(3,6);
t477 = -t204 / 0.4e1;
t205 = rSges(3,1) * t471 - Icges(3,5);
t476 = t205 / 0.4e1;
t475 = -t219 / 0.2e1;
t452 = rSges(3,2) * t272;
t197 = rSges(3,1) * t266 + t452;
t145 = -t197 * t375 - t253 * t494;
t474 = m(3) * t145;
t451 = rSges(3,2) * t274;
t198 = rSges(3,1) * t268 + t451;
t146 = -t198 * t373 - t253 * t495;
t473 = m(3) * t146;
t450 = rSges(3,2) * t276;
t199 = rSges(3,1) * t270 + t450;
t147 = -t199 * t371 - t253 * t496;
t472 = m(3) * t147;
t470 = pkin(2) * t266;
t469 = pkin(2) * t268;
t468 = pkin(2) * t270;
t290 = pkin(2) ^ 2;
t220 = t284 ^ 2 + t290;
t288 = pkin(3) ^ 2;
t428 = t266 * t64;
t344 = pkin(3) * t428;
t351 = pkin(3) * t484;
t461 = (-t284 * t344 + (t247 * t288 + t272 * t351 + t220) * t58) * t58;
t427 = t268 * t65;
t343 = pkin(3) * t427;
t460 = (-t284 * t343 + (t248 * t288 + t274 * t351 + t220) * t59) * t59;
t426 = t270 * t66;
t342 = pkin(3) * t426;
t459 = (-t284 * t342 + (t249 * t288 + t276 * t351 + t220) * t60) * t60;
t103 = t299 * t235 + t238 * t375;
t133 = -t163 * t251 * t235 + t238 * t253;
t398 = t160 * t251;
t188 = pkin(2) * t273 + t356;
t376 = t251 * t266;
t296 = pkin(3) * t376 - t185 * t253;
t127 = -t188 * t252 - t296 * t250;
t130 = t188 * t250 - t296 * t252;
t88 = t151 * t238 + (t127 * t229 + t130 * t223) * t235;
t94 = -t127 * t223 + t130 * t229;
t67 = -(t103 * t226 - t124 * t232) * t467 + (-t226 * t88 + t232 * t94) * t272 - (t133 * t226 + t232 * t398) * t470;
t449 = t121 * t67;
t68 = (t103 * t232 + t124 * t226) * t467 + (t226 * t94 + t232 * t88) * t272 + (t133 * t232 - t226 * t398) * t470;
t448 = t121 * t68;
t104 = t298 * t236 + t239 * t373;
t134 = -t164 * t251 * t236 + t239 * t253;
t397 = t161 * t251;
t189 = pkin(2) * t275 + t354;
t374 = t251 * t268;
t295 = pkin(3) * t374 - t186 * t253;
t128 = -t189 * t252 - t295 * t250;
t131 = t189 * t250 - t295 * t252;
t89 = t152 * t239 + (t128 * t230 + t131 * t224) * t236;
t95 = -t128 * t224 + t131 * t230;
t69 = -(t104 * t227 - t125 * t233) * t466 + (-t227 * t89 + t233 * t95) * t274 - (t134 * t227 + t233 * t397) * t469;
t447 = t122 * t69;
t70 = (t104 * t233 + t125 * t227) * t466 + (t227 * t95 + t233 * t89) * t274 + (t134 * t233 - t227 * t397) * t469;
t446 = t122 * t70;
t105 = t297 * t237 + t240 * t371;
t135 = -t165 * t251 * t237 + t240 * t253;
t396 = t162 * t251;
t190 = pkin(2) * t277 + t352;
t372 = t251 * t270;
t294 = pkin(3) * t372 - t187 * t253;
t129 = -t190 * t252 - t294 * t250;
t132 = t190 * t250 - t294 * t252;
t90 = t153 * t240 + (t129 * t231 + t132 * t225) * t237;
t96 = -t129 * t225 + t132 * t231;
t71 = -(t105 * t228 - t126 * t234) * t465 + (-t228 * t90 + t234 * t96) * t276 - (t135 * t228 + t234 * t396) * t468;
t445 = t123 * t71;
t72 = (t105 * t234 + t126 * t228) * t465 + (t228 * t96 + t234 * t90) * t276 + (t135 * t234 - t228 * t396) * t468;
t444 = t123 * t72;
t443 = t136 * t76;
t442 = t136 * t79;
t441 = t137 * t77;
t440 = t137 * t80;
t439 = t138 * t78;
t438 = t138 * t81;
t437 = t139 * t82;
t436 = t139 * t85;
t435 = t140 * t83;
t434 = t140 * t86;
t433 = t141 * t84;
t432 = t141 * t87;
t431 = t238 * t97;
t430 = t239 * t98;
t429 = t240 * t99;
t425 = t273 * t58;
t424 = t275 * t59;
t423 = t277 * t60;
t422 = t284 * t58;
t421 = t284 * t59;
t420 = t284 * t60;
t203 = t278 ^ 2 + t286 + t290;
t242 = t482 * t484;
t312 = (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,1) + Icges(2,3);
t346 = -0.2e1 * rSges(3,2) * pkin(2);
t112 = t200 * t247 + (t266 * t483 + t242) * t272 + (t266 * t346 + t203) * m(3) + t312;
t419 = t112 * t136;
t113 = t200 * t248 + (t268 * t483 + t242) * t274 + (t268 * t346 + t203) * m(3) + t312;
t418 = t113 * t137;
t114 = t200 * t249 + (t270 * t483 + t242) * t276 + (t270 * t346 + t203) * m(3) + t312;
t417 = t114 * t138;
t246 = m(1) + m(2) + m(3);
t416 = t121 * t246;
t415 = t122 * t246;
t414 = t123 * t246;
t154 = -t204 * t272 - t266 * t205;
t413 = t136 * t154;
t155 = -t204 * t274 - t268 * t205;
t412 = t137 * t155;
t156 = -t204 * t276 - t270 * t205;
t411 = t138 * t156;
t212 = m(2) * rSges(2,1) + t285;
t166 = -m(3) * t494 + t212;
t202 = m(2) * rSges(2,2) - t471;
t142 = t166 * t273 - t202 * t267;
t407 = t142 * t251;
t167 = -m(3) * t495 + t212;
t143 = t167 * t275 - t202 * t269;
t406 = t143 * t251;
t168 = -m(3) * t496 + t212;
t144 = t168 * t277 - t202 * t271;
t405 = t144 * t251;
t201 = (t286 + t287) * m(3) + Icges(3,3);
t392 = t201 * t289;
t382 = t238 * t251;
t380 = t239 * t251;
t378 = t240 * t251;
t369 = t251 * t273;
t367 = t251 * t275;
t365 = t251 * t277;
t358 = t253 * t289;
t347 = -t345 / 0.2e1 + Icges(3,4) / 0.2e1;
t340 = t121 * t474;
t339 = t122 * t473;
t338 = t123 * t472;
t337 = t266 * t422;
t336 = t268 * t421;
t335 = t270 * t420;
t334 = t100 * t139 * t238;
t333 = t101 * t140 * t239;
t332 = t102 * t141 * t240;
t331 = t121 * t407;
t330 = t122 * t406;
t329 = t123 * t405;
t328 = t136 * t407;
t327 = t137 * t406;
t326 = t138 * t405;
t325 = t154 * t410;
t324 = t139 * t392;
t323 = t155 * t409;
t322 = t140 * t392;
t321 = t156 * t408;
t320 = t141 * t392;
t222 = rSges(3,2) * t341;
t318 = t410 * t474;
t317 = t409 * t473;
t316 = t408 * t472;
t315 = t289 * t334;
t314 = t289 * t333;
t313 = t289 * t332;
t304 = t106 * t250 + t107 * t252;
t302 = t108 * t250 + t109 * t252;
t300 = t110 * t250 + t111 * t252;
t293 = t488 * t267 + t304 * t273;
t292 = t489 * t269 + t302 * t275;
t291 = t490 * t271 + t300 * t277;
t265 = xDDP(1);
t264 = xDDP(2);
t263 = xDDP(3);
t184 = (t287 / 0.2e1 - t286 / 0.2e1) * m(3) - Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1;
t75 = -((-t162 * t277 + t165 * t359) * t240 - t237 * t371) * t465 + ((t162 * t190 + t294 * t165) * t240 + t153 * t237) * t276 + (t165 * t378 + t237 * t253) * t468;
t74 = -((-t161 * t275 + t164 * t361) * t239 - t236 * t373) * t466 + ((t161 * t189 + t295 * t164) * t239 + t152 * t236) * t274 + (t164 * t380 + t236 * t253) * t469;
t73 = -((-t160 * t273 + t163 * t363) * t238 - t235 * t375) * t467 + ((t160 * t188 + t296 * t163) * t238 + t151 * t235) * t272 + (t163 * t382 + t235 * t253) * t470;
t63 = t66 ^ 2;
t62 = t65 ^ 2;
t61 = t64 ^ 2;
t54 = t201 * t313 + (-t156 * t429 + t75 * t472) * t123;
t53 = t201 * t314 + (-t155 * t430 + t74 * t473) * t122;
t52 = t201 * t315 + (-t154 * t431 + t73 * t474) * t121;
t51 = t313 * t472 + (-t144 * t99 * t378 + t246 * t75) * t123;
t50 = t314 * t473 + (-t143 * t98 * t380 + t246 * t74) * t122;
t49 = t315 * t474 + (-t142 * t97 * t382 + t246 * t73) * t121;
t48 = t156 * t313 + (-t114 * t429 + t75 * t405) * t123;
t47 = t155 * t314 + (-t113 * t430 + t74 * t406) * t122;
t46 = t154 * t315 + (-t112 * t431 + t73 * t407) * t121;
t45 = t87 * t320 + t72 * t338 + t78 * t411;
t44 = t84 * t320 + t71 * t338 + t81 * t411;
t43 = t86 * t322 + t70 * t339 + t77 * t412;
t42 = t83 * t322 + t69 * t339 + t80 * t412;
t41 = t85 * t324 + t68 * t340 + t76 * t413;
t40 = t82 * t324 + t67 * t340 + t79 * t413;
t39 = t87 * t316 + t78 * t326 + t72 * t414;
t38 = t84 * t316 + t81 * t326 + t71 * t414;
t37 = t86 * t317 + t77 * t327 + t70 * t415;
t36 = t83 * t317 + t80 * t327 + t69 * t415;
t35 = t85 * t318 + t76 * t328 + t68 * t416;
t34 = t82 * t318 + t79 * t328 + t67 * t416;
t33 = t87 * t321 + t72 * t329 + t78 * t417;
t32 = t84 * t321 + t71 * t329 + t81 * t417;
t31 = t86 * t323 + t70 * t330 + t77 * t418;
t30 = t83 * t323 + t69 * t330 + t80 * t418;
t29 = t85 * t325 + t68 * t331 + t76 * t419;
t28 = t82 * t325 + t67 * t331 + t79 * t419;
t24 = t335 - t479;
t23 = t336 - t480;
t22 = t337 - t481;
t18 = (-t276 * t459 - (pkin(2) * t66 - t24 * t276) * t479) * t123;
t17 = (-t274 * t460 - (pkin(2) * t65 - t23 * t274) * t480) * t122;
t16 = (-t272 * t461 - (pkin(2) * t64 - t22 * t272) * t481) * t121;
t15 = t123 * t358 * t459 + (-t253 * t335 + (-t159 * t372 + (pkin(2) * t276 + t465) * t253) * t66) * t138 * t66;
t14 = t122 * t358 * t460 + (-t253 * t336 + (-t158 * t374 + (pkin(2) * t274 + t466) * t253) * t65) * t137 * t65;
t13 = t121 * t358 * t461 + (-t253 * t337 + (-t157 * t376 + (pkin(2) * t272 + t467) * t253) * t64) * t136 * t64;
t12 = (((t253 * t66 + t60 * t365) * t465 + ((-t342 + t420) * t271 + pkin(2) * t423) * t366 + t253 * t24) * t60 + (t66 * t365 + (t249 * t253 - t353 * t372 - t253) * t60) * t479) * t123;
t11 = (((t253 * t65 + t59 * t367) * t466 + ((-t343 + t421) * t269 + pkin(2) * t424) * t368 + t253 * t23) * t59 + (t65 * t367 + (t248 * t253 - t355 * t374 - t253) * t59) * t480) * t122;
t10 = (((t253 * t64 + t58 * t369) * t467 + ((-t344 + t422) * t267 + pkin(2) * t425) * t370 + t253 * t22) * t58 + (t64 * t369 + (t247 * t253 - t357 * t376 - t253) * t58) * t481) * t121;
t9 = -t18 * t472 - t156 * t12 - t201 * t15 + 0.2e1 * ((t184 * t270 + t222) * t276 + t347 + t485) * t57 + ((t493 * rSges(3,1) + t291 * rSges(3,2)) * t276 + t270 * (t291 * rSges(3,1) - t493 * rSges(3,2))) * m(3);
t8 = -t17 * t473 - t155 * t11 - t201 * t14 + 0.2e1 * ((t184 * t268 + t222) * t274 + t347 + t486) * t56 + ((t492 * rSges(3,1) + t292 * rSges(3,2)) * t274 + t268 * (t292 * rSges(3,1) - t492 * rSges(3,2))) * m(3);
t7 = -t16 * t474 - t154 * t10 - t201 * t13 + 0.2e1 * ((t184 * t266 + t222) * t272 + t347 + t487) * t55 + ((t491 * rSges(3,1) + t293 * rSges(3,2)) * t272 + t266 * (t293 * rSges(3,1) - t491 * rSges(3,2))) * m(3);
t6 = -t18 * t405 - t114 * t12 - t156 * t15 - 0.4e1 * ((t477 * t270 + t476 * t276) * t66 + ((t270 * t478 + t222) * t276 + t475 + t485) * t60) * t66 + (-t168 * t490 + t300 * t202) * t277 - t271 * (-t300 * t168 - t202 * t490);
t5 = -t17 * t406 - t113 * t11 - t155 * t14 - 0.4e1 * ((t477 * t268 + t476 * t274) * t65 + ((t268 * t478 + t222) * t274 + t475 + t486) * t59) * t65 + (-t167 * t489 + t302 * t202) * t275 - t269 * (-t302 * t167 - t202 * t489);
t4 = -t16 * t407 - t112 * t10 - t154 * t13 - 0.4e1 * ((t477 * t266 + t476 * t272) * t64 + ((t266 * t478 + t222) * t272 + t475 + t487) * t58) * t64 + (-t166 * t488 + t304 * t202) * t273 - t267 * (-t304 * t166 - t202 * t488);
t3 = (-t144 * t12 + (-t202 * t277 - t212 * t271) * t57) * t251 + (-t18 - t150) * t246 + (-t147 * t15 + (-0.2e1 * (rSges(3,1) * t426 + t66 * t450) * t423 + t496 * t271 * (t57 + t63)) * t251 - t63 * t253 * t199) * m(3);
t2 = (-t143 * t11 + (-t202 * t275 - t212 * t269) * t56) * t251 + (-t17 - t149) * t246 + (-t146 * t14 + (-0.2e1 * (rSges(3,1) * t427 + t65 * t451) * t424 + t495 * t269 * (t56 + t62)) * t251 - t62 * t253 * t198) * m(3);
t1 = (-t10 * t142 + (-t202 * t273 - t212 * t267) * t55) * t251 + (-t16 - t148) * t246 + (-t145 * t13 + (-0.2e1 * (rSges(3,1) * t428 + t64 * t452) * t425 + t494 * t267 * (t55 + t61)) * t251 - t61 * t253 * t197) * m(3);
t19 = [-m(4) * g(1) + (t75 * t3 - t6 * t429) * t123 + (t74 * t2 - t5 * t430) * t122 + (t73 * t1 - t4 * t431) * t121 + (t9 * t332 + t8 * t333 + t7 * t334) * t289 + (t49 * t449 + t50 * t447 + t51 * t445 + t46 * t442 + t47 * t440 + t48 * t438 + (t54 * t433 + t53 * t435 + t52 * t437) * t289) * t264 + (t49 * t448 + t50 * t446 + t51 * t444 + t46 * t443 + t47 * t441 + t48 * t439 + (t54 * t432 + t53 * t434 + t52 * t436) * t289) * t263 + (m(4) + (-t48 * t429 + t51 * t75) * t123 + (-t47 * t430 + t50 * t74) * t122 + (-t46 * t431 + t49 * t73) * t121 + (t54 * t332 + t53 * t333 + t52 * t334) * t289) * t265; t1 * t449 + t2 * t447 + t3 * t445 + t4 * t442 + t5 * t440 + t6 * t438 - m(4) * g(2) + (t9 * t433 + t8 * t435 + t7 * t437) * t289 + ((-t32 * t429 + t38 * t75) * t123 + (-t30 * t430 + t36 * t74) * t122 + (-t28 * t431 + t34 * t73) * t121 + (t44 * t332 + t42 * t333 + t40 * t334) * t289) * t265 + (t34 * t448 + t36 * t446 + t38 * t444 + t28 * t443 + t30 * t441 + t32 * t439 + (t40 * t436 + t42 * t434 + t44 * t432) * t289) * t263 + (t34 * t449 + t36 * t447 + t38 * t445 + t28 * t442 + t30 * t440 + t32 * t438 + m(4) + (t40 * t437 + t42 * t435 + t44 * t433) * t289) * t264; t1 * t448 + t2 * t446 + t3 * t444 + t4 * t443 + t5 * t441 + t6 * t439 - m(4) * g(3) + (t9 * t432 + t8 * t434 + t7 * t436) * t289 + ((-t33 * t429 + t39 * t75) * t123 + (-t31 * t430 + t37 * t74) * t122 + (-t29 * t431 + t35 * t73) * t121 + (t45 * t332 + t43 * t333 + t41 * t334) * t289) * t265 + (t35 * t449 + t37 * t447 + t39 * t445 + t29 * t442 + t31 * t440 + t33 * t438 + (t41 * t437 + t43 * t435 + t45 * t433) * t289) * t264 + (t35 * t448 + t37 * t446 + t39 * t444 + t29 * t443 + t31 * t441 + t33 * t439 + m(4) + (t41 * t436 + t43 * t434 + t45 * t432) * t289) * t263;];
tauX  = t19;
