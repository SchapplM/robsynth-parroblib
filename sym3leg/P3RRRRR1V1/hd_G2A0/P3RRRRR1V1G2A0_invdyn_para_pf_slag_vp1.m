% Calculate vector of inverse dynamics forces for parallel robot
% P3RRRRR1V1G2A0
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1]';
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
% Datum: 2020-08-07 03:36
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRRRR1V1G2A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR1V1G2A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRRRR1V1G2A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRRRR1V1G2A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR1V1G2A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR1V1G2A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRRRR1V1G2A0_invdyn_para_pf_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR1V1G2A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRRRR1V1G2A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRRRR1V1G2A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR1V1G2A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR1V1G2A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:33:20
% EndTime: 2020-08-07 03:33:31
% DurationCPUTime: 11.28s
% Computational Cost: add. (42192->677), mult. (75990->1110), div. (16929->8), fcn. (76497->54), ass. (0->456)
t242 = qJ(2,3) + qJ(3,3);
t133 = cos(qJ(1,3) - t242) + cos(qJ(1,3) + t242);
t244 = qJ(2,2) + qJ(3,2);
t134 = cos(qJ(1,2) - t244) + cos(qJ(1,2) + t244);
t246 = qJ(2,1) + qJ(3,1);
t135 = cos(qJ(1,1) - t246) + cos(qJ(1,1) + t246);
t262 = cos(qJ(3,3));
t232 = t262 ^ 2;
t265 = cos(qJ(3,2));
t235 = t265 ^ 2;
t268 = cos(qJ(3,1));
t238 = t268 ^ 2;
t264 = cos(qJ(1,3));
t518 = 0.2e1 * t264 ^ 2;
t267 = cos(qJ(1,2));
t517 = 0.2e1 * t267 ^ 2;
t270 = cos(qJ(1,1));
t516 = 0.2e1 * t270 ^ 2;
t278 = xDP(3);
t259 = sin(qJ(3,1));
t231 = 0.1e1 / t259;
t292 = 0.1e1 / pkin(3);
t294 = 0.1e1 / pkin(2);
t393 = t292 * t294;
t336 = t231 * t393;
t467 = pkin(3) * t268;
t187 = pkin(2) + t467;
t269 = cos(qJ(2,1));
t156 = t187 * t269;
t260 = sin(qJ(2,1));
t397 = t259 * t260;
t361 = pkin(3) * t397;
t123 = t156 - t361;
t507 = t123 * t270;
t312 = t336 * t507;
t300 = t278 * t312;
t280 = xDP(1);
t396 = t259 * t269;
t360 = pkin(3) * t396;
t120 = t187 * t260 + t360;
t249 = legFrame(1,2);
t213 = sin(t249);
t216 = cos(t249);
t261 = sin(qJ(1,1));
t409 = t216 * t261;
t99 = t120 * t213 - t123 * t409;
t325 = t99 * t336;
t78 = t280 * t325;
t279 = xDP(2);
t412 = t213 * t261;
t102 = t120 * t216 + t123 * t412;
t322 = t102 * t336;
t81 = t279 * t322;
t63 = t78 + t81 - t300;
t406 = t231 * t294;
t332 = t406 / 0.2e1;
t319 = t135 * t332;
t303 = -t268 * t269 + t397;
t304 = -t260 * t268 - t396;
t96 = t304 * t213 - t303 * t409;
t351 = t96 * t406;
t95 = t304 * t216 + t303 * t412;
t352 = t95 * t406;
t66 = t278 * t319 + t279 * t352 + t280 * t351;
t42 = t63 + t66;
t488 = pkin(3) * t42;
t515 = -0.2e1 * t488;
t256 = sin(qJ(3,2));
t230 = 0.1e1 / t256;
t337 = t230 * t393;
t468 = pkin(3) * t265;
t186 = pkin(2) + t468;
t266 = cos(qJ(2,2));
t155 = t186 * t266;
t257 = sin(qJ(2,2));
t401 = t256 * t257;
t363 = pkin(3) * t401;
t122 = t155 - t363;
t508 = t122 * t267;
t313 = t337 * t508;
t301 = t278 * t313;
t400 = t256 * t266;
t362 = pkin(3) * t400;
t119 = t186 * t257 + t362;
t248 = legFrame(2,2);
t212 = sin(t248);
t215 = cos(t248);
t258 = sin(qJ(1,2));
t410 = t215 * t258;
t98 = t119 * t212 - t122 * t410;
t326 = t98 * t337;
t77 = t280 * t326;
t413 = t212 * t258;
t101 = t119 * t215 + t122 * t413;
t323 = t101 * t337;
t80 = t279 * t323;
t62 = t77 + t80 - t301;
t407 = t230 * t294;
t333 = t407 / 0.2e1;
t320 = t134 * t333;
t305 = -t265 * t266 + t401;
t306 = -t257 * t265 - t400;
t94 = t306 * t212 - t305 * t410;
t353 = t94 * t407;
t93 = t306 * t215 + t305 * t413;
t354 = t93 * t407;
t65 = t278 * t320 + t279 * t354 + t280 * t353;
t41 = t62 + t65;
t489 = pkin(3) * t41;
t514 = -0.2e1 * t489;
t253 = sin(qJ(3,3));
t229 = 0.1e1 / t253;
t338 = t229 * t393;
t469 = pkin(3) * t262;
t185 = pkin(2) + t469;
t263 = cos(qJ(2,3));
t154 = t185 * t263;
t254 = sin(qJ(2,3));
t405 = t253 * t254;
t365 = pkin(3) * t405;
t121 = t154 - t365;
t509 = t121 * t264;
t314 = t338 * t509;
t302 = t278 * t314;
t404 = t253 * t263;
t364 = pkin(3) * t404;
t118 = t185 * t254 + t364;
t247 = legFrame(3,2);
t211 = sin(t247);
t214 = cos(t247);
t255 = sin(qJ(1,3));
t411 = t214 * t255;
t97 = t118 * t211 - t121 * t411;
t327 = t97 * t338;
t76 = t280 * t327;
t414 = t211 * t255;
t100 = t118 * t214 + t121 * t414;
t324 = t100 * t338;
t79 = t279 * t324;
t61 = t76 + t79 - t302;
t408 = t229 * t294;
t334 = t408 / 0.2e1;
t321 = t133 * t334;
t307 = -t262 * t263 + t405;
t308 = -t254 * t262 - t404;
t92 = t308 * t211 - t307 * t411;
t355 = t92 * t408;
t91 = t308 * t214 + t307 * t414;
t356 = t91 * t408;
t64 = t278 * t321 + t279 * t356 + t280 * t355;
t40 = t61 + t64;
t490 = pkin(3) * t40;
t513 = -0.2e1 * t490;
t494 = m(3) * rSges(3,3);
t149 = rSges(3,1) * t253 + rSges(3,2) * t262;
t512 = t149 * t61;
t150 = rSges(3,1) * t256 + rSges(3,2) * t265;
t511 = t150 * t62;
t151 = rSges(3,1) * t259 + rSges(3,2) * t268;
t510 = t151 * t63;
t286 = 0.2e1 * qJ(2,1);
t245 = t286 + qJ(3,1);
t203 = sin(t245);
t209 = cos(t245);
t506 = rSges(3,1) * t203 + rSges(3,2) * t209;
t285 = 0.2e1 * qJ(2,2);
t243 = qJ(3,2) + t285;
t201 = sin(t243);
t207 = cos(t243);
t505 = rSges(3,1) * t201 + rSges(3,2) * t207;
t284 = 0.2e1 * qJ(2,3);
t241 = t284 + qJ(3,3);
t199 = sin(t241);
t205 = cos(t241);
t504 = rSges(3,1) * t199 + rSges(3,2) * t205;
t498 = m(3) / 0.2e1;
t503 = 0.2e1 * pkin(2) * t498;
t502 = -0.2e1 * t64;
t501 = -0.2e1 * t65;
t500 = -0.2e1 * t66;
t291 = pkin(3) ^ 2;
t499 = -0.2e1 * t291;
t497 = m(3) * pkin(1);
t281 = pkin(2) * m(3);
t276 = m(2) * rSges(2,2);
t496 = m(2) * rSges(2,3);
t495 = m(3) * rSges(3,2);
t493 = pkin(2) * t64;
t492 = pkin(2) * t65;
t491 = pkin(2) * t66;
t487 = t133 / 0.2e1;
t486 = t134 / 0.2e1;
t485 = t135 / 0.2e1;
t288 = rSges(2,2) ^ 2;
t290 = rSges(2,1) ^ 2;
t293 = pkin(2) ^ 2;
t143 = t293 * m(3) + (-t288 + t290) * m(2) + Icges(2,2) - Icges(2,1);
t484 = t143 / 0.2e1;
t287 = rSges(3,2) ^ 2;
t289 = rSges(3,1) ^ 2;
t157 = m(3) * (-t287 + t289) - Icges(3,1) + Icges(3,2);
t483 = t157 / 0.2e1;
t183 = m(2) * rSges(2,1) + t281;
t482 = pkin(1) * t183;
t200 = sin(t242);
t481 = pkin(1) * t200;
t202 = sin(t244);
t480 = pkin(1) * t202;
t204 = sin(t246);
t479 = pkin(1) * t204;
t206 = cos(t242);
t478 = pkin(1) * t206;
t208 = cos(t244);
t477 = pkin(1) * t208;
t210 = cos(t246);
t476 = pkin(1) * t210;
t475 = pkin(2) * t199;
t474 = pkin(2) * t201;
t473 = pkin(2) * t203;
t472 = pkin(3) * t232;
t471 = pkin(3) * t235;
t470 = pkin(3) * t238;
t220 = t262 * pkin(2);
t221 = t265 * pkin(2);
t222 = t268 * pkin(2);
t466 = t40 * t61;
t464 = t41 * t62;
t462 = t42 * t63;
t451 = rSges(3,1) * t262;
t450 = rSges(3,1) * t265;
t449 = rSges(3,1) * t268;
t445 = t292 * t97;
t444 = t292 * t98;
t443 = t292 * t99;
t442 = t100 * t292;
t441 = t101 * t292;
t440 = t102 * t292;
t191 = rSges(3,2) * t494 - Icges(3,6);
t192 = rSges(3,1) * t494 - Icges(3,5);
t109 = t191 * t206 + t192 * t200;
t439 = t109 * t292;
t110 = t191 * t208 + t192 * t202;
t438 = t110 * t292;
t111 = t191 * t210 + t192 * t204;
t437 = t111 * t292;
t377 = t253 * t495;
t172 = pkin(2) * t377;
t217 = t287 + t289;
t371 = pkin(2) * t451;
t115 = Icges(3,3) - t172 + (t217 + t371) * m(3);
t436 = t115 * t292;
t376 = t256 * t495;
t173 = pkin(2) * t376;
t369 = pkin(2) * t450;
t116 = Icges(3,3) - t173 + (t217 + t369) * m(3);
t435 = t116 * t292;
t375 = t259 * t495;
t174 = pkin(2) * t375;
t367 = pkin(2) * t449;
t117 = Icges(3,3) - t174 + (t217 + t367) * m(3);
t434 = t117 * t292;
t130 = 0.1e1 / (pkin(2) * t263 + pkin(3) * t206 + pkin(1));
t433 = t130 * t255;
t432 = t130 * t264;
t131 = 0.1e1 / (pkin(2) * t266 + pkin(3) * t208 + pkin(1));
t431 = t131 * t258;
t430 = t131 * t267;
t132 = 0.1e1 / (pkin(2) * t269 + pkin(3) * t210 + pkin(1));
t429 = t132 * t261;
t428 = t132 * t270;
t223 = sin(t284);
t427 = t143 * t223;
t224 = sin(t285);
t426 = t143 * t224;
t225 = sin(t286);
t425 = t143 * t225;
t188 = 0.2e1 * t242;
t165 = sin(t188);
t424 = t157 * t165;
t189 = 0.2e1 * t244;
t166 = sin(t189);
t423 = t157 * t166;
t190 = 0.2e1 * t246;
t167 = sin(t190);
t422 = t157 * t167;
t162 = m(3) * t217 + Icges(3,3);
t421 = t162 * t292;
t168 = cos(t188);
t193 = rSges(3,1) * t495 - Icges(3,4);
t420 = t193 * t168;
t169 = cos(t189);
t419 = t193 * t169;
t170 = cos(t190);
t418 = t193 * t170;
t195 = rSges(2,1) * t276 - Icges(2,4);
t226 = cos(t284);
t417 = t195 * t226;
t227 = cos(t285);
t416 = t195 * t227;
t228 = cos(t286);
t415 = t195 * t228;
t403 = t254 * t263;
t402 = t255 * t264;
t399 = t257 * t266;
t398 = t258 * t267;
t395 = t260 * t269;
t394 = t261 * t270;
t219 = t288 + t290;
t392 = t291 - t293;
t391 = 0.2e1 * pkin(1);
t390 = pkin(1) * t276;
t386 = pkin(2) * t494;
t385 = rSges(2,1) * t496;
t384 = rSges(2,2) * t496;
t383 = pkin(3) * t220;
t382 = pkin(3) * t221;
t381 = pkin(3) * t222;
t380 = t254 * t502;
t379 = t257 * t501;
t378 = t260 * t500;
t374 = t40 * t469;
t373 = t41 * t468;
t372 = t42 * t467;
t370 = t64 * t220;
t368 = t65 * t221;
t366 = t66 * t222;
t359 = t493 / 0.2e1;
t358 = t492 / 0.2e1;
t357 = t491 / 0.2e1;
t350 = 0.2e1 * rSges(3,1) * t497;
t349 = -0.2e1 * pkin(1) * t495;
t348 = 0.2e1 * t281;
t347 = t292 * t509;
t346 = t292 * t508;
t345 = t292 * t507;
t344 = t211 * t432;
t343 = t214 * t432;
t342 = t212 * t430;
t341 = t215 * t430;
t340 = t213 * t428;
t339 = t216 * t428;
t184 = m(1) * rSges(1,1) + m(2) * pkin(1);
t335 = -t184 - t497;
t175 = m(1) * rSges(1,2) + t496;
t331 = (t175 + t494) * g(3);
t318 = m(2) * t219 + Icges(2,3) + Icges(3,3);
t139 = g(1) * t214 - g(2) * t211;
t317 = g(3) * t264 + t139 * t255;
t140 = g(1) * t215 - g(2) * t212;
t316 = g(3) * t267 + t140 * t258;
t141 = g(1) * t216 - g(2) * t213;
t315 = g(3) * t270 + t141 * t261;
t124 = -m(3) * t451 - t183 + t377;
t127 = t149 * m(3) + t276;
t311 = -t124 * t263 - t127 * t254;
t125 = -m(3) * t450 - t183 + t376;
t128 = t150 * m(3) + t276;
t310 = -t266 * t125 - t128 * t257;
t126 = -m(3) * t449 - t183 + t375;
t129 = t151 * m(3) + t276;
t309 = -t269 * t126 - t129 * t260;
t299 = pkin(1) * t405 - pkin(3) + t472;
t298 = pkin(1) * t401 - pkin(3) + t471;
t297 = pkin(1) * t397 - pkin(3) + t470;
t283 = 0.2e1 * pkin(1) ^ 2;
t296 = rSges(3,3) ^ 2 + t283 / 0.2e1 + t287 / 0.2e1 + t289 / 0.2e1 + t293 / 0.2e1;
t295 = Icges(1,3) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + Icges(2,1) / 0.2e1 + Icges(3,1) / 0.2e1 + Icges(2,2) / 0.2e1 + Icges(3,2) / 0.2e1 + ((2 * rSges(2,3) ^ 2) + t219 + t283) * m(2) / 0.2e1;
t282 = pkin(1) * g(3);
t252 = xDDP(1);
t251 = xDDP(2);
t250 = xDDP(3);
t239 = t269 ^ 2;
t236 = t266 ^ 2;
t233 = t263 ^ 2;
t194 = -Icges(2,6) + t384;
t182 = t293 + t217;
t171 = t384 / 0.2e1 - Icges(2,6) / 0.2e1;
t164 = t184 * g(3);
t161 = pkin(1) - 0.2e1 * t361;
t160 = pkin(1) - 0.2e1 * t363;
t159 = pkin(1) - 0.2e1 * t365;
t152 = -Icges(2,5) + t385 + t386;
t148 = -pkin(3) + t222 + 0.2e1 * t470;
t147 = -pkin(3) + t221 + 0.2e1 * t471;
t146 = -pkin(3) + t220 + 0.2e1 * t472;
t144 = -t385 / 0.2e1 - t386 / 0.2e1 + Icges(2,5) / 0.2e1;
t138 = g(1) * t213 + g(2) * t216;
t137 = g(1) * t212 + g(2) * t215;
t136 = g(1) * t211 + g(2) * t214;
t108 = -0.2e1 * t174 + (t182 + 0.2e1 * t367) * m(3) + t318;
t107 = -0.2e1 * t173 + (t182 + 0.2e1 * t369) * m(3) + t318;
t106 = -0.2e1 * t172 + (t182 + 0.2e1 * t371) * m(3) + t318;
t105 = t152 * t260 + t194 * t269 + t111;
t104 = t152 * t257 + t194 * t266 + t110;
t103 = t152 * t254 + t194 * t263 + t109;
t75 = (-t261 * t278 + (-t213 * t279 + t216 * t280) * t270) * t132;
t74 = (-t258 * t278 + (-t212 * t279 + t215 * t280) * t267) * t131;
t73 = (-t255 * t278 + (-t211 * t279 + t214 * t280) * t264) * t130;
t72 = t75 ^ 2;
t71 = t74 ^ 2;
t70 = t73 ^ 2;
t69 = -t111 * t429 + (t117 * t485 - t162 * t345) * t406;
t68 = -t110 * t431 + (t116 * t486 - t162 * t346) * t407;
t67 = -t109 * t433 + (t115 * t487 - t162 * t347) * t408;
t60 = -t105 * t429 + (t108 * t485 - t117 * t345) * t406;
t59 = -t104 * t431 + (t107 * t486 - t116 * t346) * t407;
t58 = -t103 * t433 + (t106 * t487 - t115 * t347) * t408;
t57 = t170 * t483 + t228 * t484 - t174 + (t183 * t269 - t260 * t276) * t391 - t193 * t167 - t195 * t225 + ((-t473 - 0.2e1 * t479) * rSges(3,2) + (0.2e1 * t476 + (t209 + t268) * pkin(2)) * rSges(3,1) + t296) * m(3) + t295;
t56 = t169 * t483 + t227 * t484 - t173 + (t183 * t266 - t257 * t276) * t391 - t193 * t166 - t195 * t224 + ((-t474 - 0.2e1 * t480) * rSges(3,2) + (0.2e1 * t477 + (t207 + t265) * pkin(2)) * rSges(3,1) + t296) * m(3) + t295;
t55 = t168 * t483 + t226 * t484 - t172 + (t183 * t263 - t254 * t276) * t391 - t193 * t165 - t195 * t223 + ((-t475 - 0.2e1 * t481) * rSges(3,2) + (0.2e1 * t478 + (t205 + t262) * pkin(2)) * rSges(3,1) + t296) * m(3) + t295;
t54 = t111 * t339 + (t117 * t96 + t99 * t421) * t406;
t53 = t110 * t341 + (t116 * t94 + t98 * t421) * t407;
t52 = t109 * t343 + (t115 * t92 + t97 * t421) * t408;
t51 = -t111 * t340 + (t102 * t421 + t117 * t95) * t406;
t50 = -t110 * t342 + (t101 * t421 + t116 * t93) * t407;
t49 = -t109 * t344 + (t100 * t421 + t115 * t91) * t408;
t48 = t105 * t339 + (t108 * t96 + t99 * t434) * t406;
t47 = t104 * t341 + (t107 * t94 + t98 * t435) * t407;
t46 = t103 * t343 + (t106 * t92 + t97 * t436) * t408;
t45 = -t105 * t340 + (t102 * t434 + t108 * t95) * t406;
t44 = -t104 * t342 + (t101 * t435 + t107 * t93) * t407;
t43 = -t103 * t344 + (t100 * t436 + t106 * t91) * t408;
t39 = t78 / 0.2e1 + t81 / 0.2e1 - t300 / 0.2e1 + t66;
t38 = t77 / 0.2e1 + t80 / 0.2e1 - t301 / 0.2e1 + t65;
t37 = t76 / 0.2e1 + t79 / 0.2e1 - t302 / 0.2e1 + t64;
t36 = -t57 * t429 + (t105 * t485 - t111 * t345) * t406;
t35 = -t56 * t431 + (t104 * t486 - t110 * t346) * t407;
t34 = -t55 * t433 + (t103 * t487 - t109 * t347) * t408;
t33 = t57 * t339 + (t105 * t96 + t99 * t437) * t406;
t32 = t56 * t341 + (t104 * t94 + t98 * t438) * t407;
t31 = t55 * t343 + (t103 * t92 + t97 * t439) * t408;
t30 = -t57 * t340 + (t102 * t437 + t105 * t95) * t406;
t29 = -t56 * t342 + (t101 * t438 + t104 * t93) * t407;
t28 = -t55 * t344 + (t100 * t439 + t103 * t91) * t408;
t27 = t372 + t491;
t26 = t373 + t492;
t25 = t374 + t493;
t18 = (pkin(2) * t378 + t204 * t515) * t75 * t132;
t17 = (pkin(2) * t379 + t202 * t514) * t74 * t131;
t16 = (pkin(2) * t380 + t200 * t513) * t73 * t130;
t15 = ((t222 + pkin(3)) * t462 + (-t72 * ((t238 * t499 - 0.2e1 * t381 + t392) * t239 - t161 * t156 + pkin(3) * t297) + (t291 * t42 + t293 * t66 + 0.2e1 * t39 * t381) * t66) * t292) * t406;
t14 = ((t221 + pkin(3)) * t464 + (-t71 * ((t235 * t499 - 0.2e1 * t382 + t392) * t236 - t160 * t155 + pkin(3) * t298) + (t291 * t41 + t293 * t65 + 0.2e1 * t38 * t382) * t65) * t292) * t407;
t13 = ((t220 + pkin(3)) * t466 + (-t70 * ((t232 * t499 - 0.2e1 * t383 + t392) * t233 - t159 * t154 + pkin(3) * t299) + (t291 * t40 + t293 * t64 + 0.2e1 * t37 * t383) * t64) * t292) * t408;
t12 = -pkin(3) * t406 * t462 + ((-0.4e1 * ((t268 * t357 + (t238 - 0.1e1 / 0.2e1) * t488) * t395 + ((t357 + t372) * t239 - t27 / 0.2e1) * t259) * t394 + (t516 - 0.2e1) * t75 * (t148 * t239 + (-pkin(2) * t397 + t161 * t268) * t269 - t297) + t135 * ((t260 * t27 + t42 * t360) * t261 - t270 * (pkin(1) + t123) * t75)) * t75 + (((t366 + (0.2e1 * t238 - 0.1e1) * t488) * t239 - (0.2e1 * t372 + t491) * t259 * t395 + t488 - t488 * t238) * t516 - 0.2e1 * t75 * (t148 * t395 + ((pkin(2) + 0.2e1 * t467) * t239 - t187) * t259) * t394 - 0.2e1 * t366 + t515 + t135 * ((-t269 * t27 + t42 * t361) * t270 + t261 * t120 * t75)) * t66) * t332;
t11 = -pkin(3) * t407 * t464 + ((-0.4e1 * ((t265 * t358 + (t235 - 0.1e1 / 0.2e1) * t489) * t399 + ((t358 + t373) * t236 - t26 / 0.2e1) * t256) * t398 + (t517 - 0.2e1) * t74 * (t147 * t236 + (-pkin(2) * t401 + t160 * t265) * t266 - t298) + t134 * ((t257 * t26 + t41 * t362) * t258 - t267 * (pkin(1) + t122) * t74)) * t74 + (((t368 + (0.2e1 * t235 - 0.1e1) * t489) * t236 - (0.2e1 * t373 + t492) * t256 * t399 + t489 - t489 * t235) * t517 - 0.2e1 * t74 * (t147 * t399 + ((pkin(2) + 0.2e1 * t468) * t236 - t186) * t256) * t398 - 0.2e1 * t368 + t514 + t134 * ((-t26 * t266 + t41 * t363) * t267 + t258 * t119 * t74)) * t65) * t333;
t10 = -pkin(3) * t408 * t466 + ((-0.4e1 * ((t262 * t359 + (t232 - 0.1e1 / 0.2e1) * t490) * t403 + ((t359 + t374) * t233 - t25 / 0.2e1) * t253) * t402 + (t518 - 0.2e1) * t73 * (t146 * t233 + (-pkin(2) * t405 + t159 * t262) * t263 - t299) + t133 * ((t25 * t254 + t40 * t364) * t255 - t264 * (pkin(1) + t121) * t73)) * t73 + (((t370 + (0.2e1 * t232 - 0.1e1) * t490) * t233 - (0.2e1 * t374 + t493) * t253 * t403 + t490 - t490 * t232) * t518 - 0.2e1 * t73 * (t146 * t403 + ((pkin(2) + 0.2e1 * t469) * t233 - t185) * t253) * t402 - 0.2e1 * t370 + t513 + t133 * ((-t25 * t263 + t40 * t365) * t264 + t255 * t118 * t73)) * t64) * t334;
t9 = -t111 * t18 - t117 * t12 - t162 * t15 + t66 ^ 2 * t151 * t503 + m(3) * ((rSges(3,1) * t138 + t315 * rSges(3,2)) * t210 + t204 * (t315 * rSges(3,1) - rSges(3,2) * t138)) + (t167 * t483 + t418 + ((rSges(3,1) * t204 + rSges(3,2) * t210) * t391 + (t151 + t506) * pkin(2)) * t498) * t72;
t8 = -t110 * t17 - t116 * t11 - t162 * t14 + t65 ^ 2 * t150 * t503 + m(3) * ((rSges(3,1) * t137 + t316 * rSges(3,2)) * t208 + t202 * (t316 * rSges(3,1) - rSges(3,2) * t137)) + (t166 * t483 + t419 + ((rSges(3,1) * t202 + rSges(3,2) * t208) * t391 + (t150 + t505) * pkin(2)) * t498) * t71;
t7 = -t109 * t16 - t115 * t10 - t162 * t13 + t64 ^ 2 * t149 * t503 + m(3) * ((rSges(3,1) * t136 + t317 * rSges(3,2)) * t206 + t200 * (t317 * rSges(3,1) - rSges(3,2) * t136)) + (t165 * t483 + t420 + ((rSges(3,1) * t200 + rSges(3,2) * t206) * t391 + (t149 + t504) * pkin(2)) * t498) * t70;
t6 = -t105 * t18 - t108 * t12 - t117 * t15 + (-t315 * t126 - t129 * t138) * t260 + t269 * (-t126 * t138 + t315 * t129) - t39 * t348 * t510 + (t422 / 0.2e1 + t418 + t425 / 0.2e1 + t415 + (t183 * t260 + t269 * t276) * pkin(1) + ((pkin(2) * t209 + t476) * rSges(3,2) + (t473 + t479) * rSges(3,1)) * m(3)) * t72;
t5 = -t104 * t17 - t107 * t11 - t116 * t14 + (-t316 * t125 - t128 * t137) * t257 + t266 * (-t125 * t137 + t316 * t128) - t38 * t348 * t511 + (t423 / 0.2e1 + t419 + t426 / 0.2e1 + t416 + (t183 * t257 + t266 * t276) * pkin(1) + ((pkin(2) * t207 + t477) * rSges(3,2) + (t474 + t480) * rSges(3,1)) * m(3)) * t71;
t4 = -t103 * t16 - t106 * t10 - t115 * t13 + (-t317 * t124 - t127 * t136) * t254 + t263 * (-t124 * t136 + t317 * t127) - t37 * t348 * t512 + (t424 / 0.2e1 + t420 + t427 / 0.2e1 + t417 + (t183 * t254 + t263 * t276) * pkin(1) + ((pkin(2) * t205 + t478) * rSges(3,2) + (t475 + t481) * rSges(3,1)) * m(3)) * t70;
t3 = -t57 * t18 - t105 * t12 - t111 * t15 + (t144 * t66 + t75 * t390) * t269 * t500 + (t171 * t66 + t75 * t482) * t378 + ((-t309 + t335) * t141 + t331) * t270 + ((rSges(3,3) * t141 + t282) * m(3) + t164 + t141 * t175 + t309 * g(3)) * t261 + (-t425 - 0.2e1 * t415) * t66 * t75 + ((-t191 * t204 + t192 * t210) * t42 + (-t350 * t204 + t349 * t210 - 0.2e1 * t418 - t422) * t75) * t42 + (-0.2e1 * t506 * t39 - t510) * t75 * t281;
t2 = -t56 * t17 - t104 * t11 - t110 * t14 + (t144 * t65 + t74 * t390) * t266 * t501 + (t171 * t65 + t74 * t482) * t379 + ((-t310 + t335) * t140 + t331) * t267 + ((rSges(3,3) * t140 + t282) * m(3) + t164 + t140 * t175 + t310 * g(3)) * t258 + (-t426 - 0.2e1 * t416) * t65 * t74 + ((-t191 * t202 + t192 * t208) * t41 + (-t350 * t202 + t349 * t208 - 0.2e1 * t419 - t423) * t74) * t41 + (-0.2e1 * t505 * t38 - t511) * t74 * t281;
t1 = -t55 * t16 - t103 * t10 - t109 * t13 + (t144 * t64 + t73 * t390) * t263 * t502 + (t171 * t64 + t73 * t482) * t380 + ((-t311 + t335) * t139 + t331) * t264 + ((rSges(3,3) * t139 + t282) * m(3) + t164 + t139 * t175 + t311 * g(3)) * t255 + (-t427 - 0.2e1 * t417) * t64 * t73 + ((-t191 * t200 + t192 * t206) * t40 + (-t350 * t200 + t349 * t206 - 0.2e1 * t420 - t424) * t73) * t40 + (-0.2e1 * t504 * t37 - t512) * t73 * t281;
t19 = [-m(4) * g(1) + t1 * t343 + t2 * t341 + t3 * t339 + t9 * t325 + t8 * t326 + t7 * t327 + t6 * t351 + t5 * t353 + t4 * t355 + (-t33 * t340 + (t54 * t440 + t48 * t95) * t406 - t32 * t342 + (t53 * t441 + t47 * t93) * t407 - t31 * t344 + (t52 * t442 + t46 * t91) * t408) * t251 + (-t33 * t429 + (-t54 * t345 + t48 * t485) * t406 - t32 * t431 + (-t53 * t346 + t47 * t486) * t407 - t31 * t433 + (-t52 * t347 + t46 * t487) * t408) * t250 + (t33 * t339 + (t54 * t443 + t48 * t96) * t406 + t32 * t341 + (t53 * t444 + t47 * t94) * t407 + t31 * t343 + (t52 * t445 + t46 * t92) * t408 + m(4)) * t252; -m(4) * g(2) - t1 * t344 - t2 * t342 - t3 * t340 + t9 * t322 + t8 * t323 + t7 * t324 + t6 * t352 + t5 * t354 + t4 * t356 + (t30 * t339 + (t51 * t443 + t45 * t96) * t406 + t29 * t341 + (t44 * t94 + t50 * t444) * t407 + t28 * t343 + (t43 * t92 + t49 * t445) * t408) * t252 + (-t30 * t429 + (-t51 * t345 + t45 * t485) * t406 - t29 * t431 + (-t50 * t346 + t44 * t486) * t407 - t28 * t433 + (-t49 * t347 + t43 * t487) * t408) * t250 + (-t30 * t340 + (t51 * t440 + t45 * t95) * t406 - t29 * t342 + (t44 * t93 + t50 * t441) * t407 - t28 * t344 + (t43 * t91 + t49 * t442) * t408 + m(4)) * t251; -m(4) * g(3) - t1 * t433 - t2 * t431 - t3 * t429 - t9 * t312 - t8 * t313 - t7 * t314 + t6 * t319 + t5 * t320 + t4 * t321 + (t36 * t339 + (t69 * t443 + t60 * t96) * t406 + t35 * t341 + (t68 * t444 + t59 * t94) * t407 + t34 * t343 + (t67 * t445 + t58 * t92) * t408) * t252 + (-t36 * t340 + (t69 * t440 + t60 * t95) * t406 - t35 * t342 + (t68 * t441 + t59 * t93) * t407 - t34 * t344 + (t67 * t442 + t58 * t91) * t408) * t251 + (-t36 * t429 + (-t69 * t345 + t60 * t485) * t406 - t35 * t431 + (-t68 * t346 + t59 * t486) * t407 - t34 * t433 + (-t67 * t347 + t58 * t487) * t408 + m(4)) * t250;];
tauX  = t19;
