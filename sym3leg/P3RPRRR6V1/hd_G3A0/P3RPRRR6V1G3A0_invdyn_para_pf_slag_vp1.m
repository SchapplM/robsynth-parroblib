% Calculate vector of inverse dynamics forces for parallel robot
% P3RPRRR6V1G3A0
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
% Datum: 2020-08-06 18:43
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRRR6V1G3A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:41:07
% EndTime: 2020-08-06 18:41:18
% DurationCPUTime: 11.08s
% Computational Cost: add. (44577->637), mult. (51672->1030), div. (7413->16), fcn. (34380->88), ass. (0->412)
t231 = xDDP(1);
t502 = t231 / 0.2e1;
t221 = sin(pkin(7));
t252 = -pkin(6) - pkin(5);
t392 = t221 * t252;
t133 = -pkin(1) + t392;
t239 = cos(qJ(1,3));
t116 = t133 * t239;
t222 = cos(pkin(7));
t233 = sin(qJ(1,3));
t390 = t233 * t252;
t395 = t221 * t233;
t501 = t116 + pkin(2) * t395 - (pkin(2) * t239 - t390) * t222;
t241 = cos(qJ(1,2));
t117 = t133 * t241;
t235 = sin(qJ(1,2));
t389 = t235 * t252;
t394 = t221 * t235;
t500 = t117 + pkin(2) * t394 - (pkin(2) * t241 - t389) * t222;
t243 = cos(qJ(1,1));
t118 = t133 * t243;
t237 = sin(qJ(1,1));
t388 = t237 * t252;
t393 = t221 * t237;
t499 = t118 + pkin(2) * t393 - (pkin(2) * t243 - t388) * t222;
t195 = qJ(1,1) + pkin(7);
t168 = cos(t195);
t229 = xDDP(3);
t228 = legFrame(1,2);
t142 = t228 + t195;
t143 = -t228 + t195;
t98 = -sin(t142) - sin(t143);
t498 = -t168 * t229 + t98 * t502;
t194 = qJ(1,2) + pkin(7);
t167 = cos(t194);
t227 = legFrame(2,2);
t140 = t227 + t194;
t141 = -t227 + t194;
t97 = -sin(t140) - sin(t141);
t497 = -t167 * t229 + t97 * t502;
t191 = qJ(1,3) + pkin(7);
t226 = legFrame(3,2);
t138 = t226 + t191;
t139 = -t226 + t191;
t96 = -sin(t138) - sin(t139);
t430 = t96 / 0.2e1;
t230 = xDDP(2);
t496 = t230 / 0.2e1;
t275 = pkin(3) ^ 2;
t276 = 0.1e1 / pkin(3);
t240 = cos(qJ(3,2));
t179 = t240 * pkin(3);
t151 = t179 + pkin(2);
t170 = t222 * pkin(1);
t120 = 0.1e1 / (t170 + t151);
t234 = sin(qJ(3,2));
t215 = 0.1e1 / t234;
t416 = t120 * t215;
t412 = t151 * t241;
t71 = (-t389 + t412) * t222 - t117 - t151 * t394;
t344 = t71 * t416;
t249 = xDP(3);
t187 = pkin(7) + qJ(3,2);
t157 = sin(t187);
t189 = -pkin(7) + qJ(3,2);
t158 = sin(t189);
t160 = sin(t194);
t447 = t235 * pkin(1);
t267 = 0.2e1 * qJ(3,2);
t201 = sin(t267);
t449 = pkin(3) * t201;
t466 = 0.2e1 * t252;
t487 = 0.2e1 * pkin(2);
t444 = (t167 * t466 + 0.2e1 * t447 + t160 * t487 + (sin(qJ(1,2) - t189) + sin(qJ(1,2) + t187)) * pkin(3)) / (t234 * t487 + t449 + (t157 + t158) * pkin(1));
t346 = t249 * t444;
t175 = cos(t227);
t251 = xDP(1);
t400 = t175 * t251;
t172 = sin(t227);
t250 = xDP(2);
t404 = t172 * t250;
t493 = (t346 / 0.4e1 - (t400 / 0.4e1 - t404 / 0.4e1) * t344) * t275 * t276;
t458 = Icges(3,2) / 0.2e1;
t492 = t458 - Icges(3,1) / 0.2e1;
t166 = cos(t191);
t99 = cos(t139) - cos(t138);
t491 = -t166 * t229 + t231 * t430 + t99 * t496;
t100 = cos(t141) - cos(t140);
t490 = t100 * t496 + t497;
t101 = cos(t143) - cos(t142);
t489 = t101 * t496 + t498;
t238 = cos(qJ(3,3));
t177 = t238 * pkin(3);
t150 = t177 + pkin(2);
t119 = 0.1e1 / (t170 + t150);
t488 = t491 * t119;
t464 = m(3) / 0.2e1;
t429 = t97 / 0.2e1;
t428 = t98 / 0.2e1;
t427 = t99 / 0.2e1;
t421 = t100 / 0.2e1;
t420 = t101 / 0.2e1;
t485 = t250 / 0.2e1;
t484 = t251 / 0.2e1;
t423 = t276 * t71;
t335 = t215 * t423;
t480 = (t400 / 0.6e1 - t404 / 0.6e1) * t335;
t479 = (t400 / 0.3e1 - t404 / 0.3e1) * t335;
t478 = (t400 / 0.2e1 - t404 / 0.2e1) * t335;
t451 = pkin(1) * t221;
t477 = -t252 + t451;
t261 = 0.2e1 * pkin(7);
t184 = cos(t261);
t278 = pkin(1) ^ 2;
t155 = t278 * t184;
t476 = -t155 - t275 / 0.2e1;
t183 = sin(t261);
t383 = t278 * t183;
t356 = 0.2e1 * t383;
t369 = -0.4e1 * t170;
t475 = t252 * t369 + t356;
t472 = 0.4e1 * pkin(2);
t471 = -0.4e1 * pkin(3);
t171 = sin(t226);
t174 = cos(t226);
t232 = sin(qJ(3,3));
t214 = 0.1e1 / t232;
t417 = t119 * t214;
t159 = sin(t191);
t264 = 0.2e1 * qJ(3,3);
t198 = sin(t264);
t365 = pkin(7) + qJ(3,3);
t367 = -pkin(7) + qJ(3,3);
t448 = t233 * pkin(1);
t445 = (t166 * t466 + 0.2e1 * t448 + t159 * t487 + (sin(qJ(1,3) - t367) + sin(qJ(1,3) + t365)) * pkin(3)) / (t232 * t487 + pkin(3) * t198 + (sin(t365) + sin(t367)) * pkin(1));
t413 = t150 * t239;
t69 = (-t390 + t413) * t222 - t116 - t150 * t395;
t39 = (t249 * t445 + (t171 * t250 - t174 * t251) * t69 * t417) * t276;
t470 = t39 ^ 2;
t40 = (t346 + (-t400 + t404) * t344) * t276;
t469 = t40 ^ 2;
t173 = sin(t228);
t176 = cos(t228);
t242 = cos(qJ(3,1));
t181 = t242 * pkin(3);
t153 = t181 + pkin(2);
t121 = 0.1e1 / (t170 + t153);
t236 = sin(qJ(3,1));
t216 = 0.1e1 / t236;
t414 = t121 * t216;
t161 = sin(t195);
t270 = 0.2e1 * qJ(3,1);
t204 = sin(t270);
t366 = pkin(7) + qJ(3,1);
t368 = -pkin(7) + qJ(3,1);
t446 = t237 * pkin(1);
t443 = (t168 * t466 + 0.2e1 * t446 + t161 * t487 + (sin(qJ(1,1) - t368) + sin(qJ(1,1) + t366)) * pkin(3)) / (t236 * t487 + pkin(3) * t204 + (sin(t366) + sin(t368)) * pkin(1));
t411 = t153 * t243;
t73 = (-t388 + t411) * t222 - t118 - t153 * t393;
t41 = (t249 * t443 + (t173 * t250 - t176 * t251) * t73 * t414) * t276;
t468 = t41 ^ 2;
t277 = pkin(2) ^ 2;
t256 = 0.2e1 * t277;
t362 = pkin(2) * t170;
t102 = 0.4e1 * t362 + t155 + t256 + t278;
t467 = -0.2e1 * t275 - 0.2e1 * t102;
t465 = -0.6e1 * t275;
t260 = -0.2e1 * t278;
t463 = pkin(1) * pkin(3);
t462 = m(2) * rSges(2,1);
t461 = m(3) * rSges(3,2);
t460 = pkin(3) * t40;
t459 = -Icges(3,5) / 0.4e1;
t457 = Icges(3,6) / 0.4e1;
t271 = rSges(3,2) ^ 2;
t272 = rSges(3,1) ^ 2;
t456 = (-t271 + t272) * t464 + t492;
t455 = t276 / 0.2e1;
t244 = rSges(3,3) + pkin(5);
t126 = rSges(3,1) * t238 - rSges(3,2) * t232;
t454 = m(3) * t126;
t127 = rSges(3,1) * t240 - rSges(3,2) * t234;
t453 = m(3) * t127;
t128 = rSges(3,1) * t242 - rSges(3,2) * t236;
t452 = m(3) * t128;
t450 = pkin(2) * t477;
t415 = t120 / 0.2e1;
t83 = t97 * t251 * t415;
t84 = t100 * t250 * t415;
t442 = t83 + t84;
t207 = cos(t264);
t316 = -0.2e1 * t252 * t170 + t383;
t89 = t316 + 0.2e1 * t450;
t441 = t207 * t89;
t213 = cos(t270);
t440 = t213 * t89;
t437 = t232 * t39;
t436 = t234 * t40;
t408 = t167 * t249;
t326 = t120 * t408;
t66 = -t326 + t442;
t435 = t234 * t66;
t434 = t236 * t41;
t68 = (t150 * t233 + t239 * t252) * t222 - t133 * t233 + t221 * t413;
t433 = t238 * t68;
t70 = (t151 * t235 + t241 * t252) * t222 - t133 * t235 + t221 * t412;
t432 = t240 * t70;
t72 = (t153 * t237 + t243 * t252) * t222 - t133 * t237 + t221 * t411;
t431 = t242 * t72;
t426 = t252 * t66;
t425 = t275 * t40;
t424 = t276 * t69;
t422 = t276 * t73;
t64 = t278 * t66;
t419 = t102 * t207;
t418 = t102 * t213;
t406 = t171 * t232;
t405 = t172 * t234;
t403 = t173 * t236;
t402 = t174 * t232;
t401 = t175 * t234;
t399 = t176 * t236;
t398 = t214 * t238;
t397 = t215 * t240;
t396 = t216 * t242;
t391 = t229 * t276;
t262 = 0.4e1 * qJ(3,3);
t387 = t275 * cos(t262);
t265 = 0.4e1 * qJ(3,2);
t386 = t275 * cos(t265);
t268 = 0.4e1 * qJ(3,1);
t385 = t275 * cos(t268);
t364 = rSges(3,1) * t461;
t382 = -t364 / 0.2e1 + Icges(3,4) / 0.2e1;
t220 = t252 ^ 2;
t381 = t220 / 0.2e1 + t278;
t266 = 0.3e1 * qJ(3,2);
t209 = cos(t266);
t380 = t209 - t240;
t379 = t220 + t278;
t378 = t271 + t272;
t377 = 0.4e1 * t460;
t353 = pkin(1) * t392;
t135 = -0.2e1 * t353;
t255 = 0.3e1 * t277;
t259 = 0.2e1 * t278;
t376 = (0.6e1 * t362 + t135 + t259 + t255 + t220 - t476) * t471;
t375 = 0.2e1 * pkin(1);
t373 = 0.2e1 * pkin(3);
t144 = t170 + pkin(2);
t371 = t144 * t471;
t370 = pkin(3) * t260;
t363 = t252 * t463;
t361 = pkin(3) * t437;
t360 = pkin(3) * t436;
t359 = pkin(3) * t434;
t358 = -0.2e1 * t477 * t275;
t357 = t144 * t465;
t149 = t177 + t487;
t355 = t149 * t170;
t152 = t181 + t487;
t354 = t152 * t170;
t217 = t238 ^ 2;
t352 = pkin(3) * (-t222 * t239 + t395) * t217;
t218 = t240 ^ 2;
t351 = pkin(3) * (-t222 * t241 + t394) * t218;
t219 = t242 ^ 2;
t350 = pkin(3) * (-t222 * t243 + t393) * t219;
t349 = t276 * t445;
t348 = t276 * t443;
t347 = t276 * t444;
t345 = g(3) * t244;
t343 = t171 * t424;
t342 = t172 * t423;
t341 = t173 * t422;
t340 = t174 * t424;
t339 = t175 * t423;
t338 = t176 * t422;
t337 = t68 * t398;
t336 = t70 * t397;
t334 = t72 * t396;
t333 = t229 * t433;
t332 = t229 * t432;
t331 = t229 * t431;
t329 = t144 * t464;
t327 = m(3) * t487;
t325 = -0.2e1 * t277 + t476;
t323 = t477 * t373;
t321 = t244 + t451;
t319 = t424 * t454;
t318 = t423 * t453;
t317 = t422 * t452;
t103 = -t321 * t461 + Icges(3,6);
t282 = -m(3) * rSges(3,1) * t321 + Icges(3,5);
t77 = t103 * t238 + t232 * t282;
t315 = t214 * t77 * t424;
t78 = t103 * t240 + t234 * t282;
t314 = t78 * t335;
t79 = t103 * t242 + t236 * t282;
t313 = t216 * t79 * t422;
t312 = t276 * t346;
t311 = rSges(3,1) * t329;
t310 = -0.2e1 * t155 - 0.4e1 * t277 + t260 - t275;
t309 = 0.3e1 / 0.8e1 * t275 + t277 / 0.2e1 + t381;
t308 = -t278 + t325;
t307 = t170 / 0.2e1 + pkin(2) / 0.2e1;
t306 = t451 / 0.4e1 + t244 / 0.4e1;
t305 = rSges(3,1) * t232 + rSges(3,2) * t238;
t304 = rSges(3,1) * t234 + rSges(3,2) * t240;
t303 = rSges(3,1) * t236 + rSges(3,2) * t242;
t65 = (-t166 * t249 + t96 * t484 + t99 * t485) * t119;
t299 = t65 * t307;
t298 = t66 * t307;
t67 = (t101 * t485 - t168 * t249 + t98 * t484) * t121;
t297 = t67 * t307;
t296 = pkin(2) + t126;
t295 = pkin(2) + t127;
t294 = pkin(2) + t128;
t290 = (rSges(2,2) * m(2) - m(3) * t244) * t221;
t263 = 0.3e1 * qJ(3,3);
t206 = cos(t263);
t289 = pkin(3) * (-pkin(2) * t238 + t144 * t206);
t269 = 0.3e1 * qJ(3,1);
t212 = cos(t269);
t288 = pkin(3) * (-pkin(2) * t242 + t144 * t212);
t283 = t64 + (0.3e1 / 0.4e1 * t83 + 0.3e1 / 0.4e1 * t84 - 0.3e1 / 0.4e1 * t326) * t275 + 0.3e1 * t66 * (t277 + t220 / 0.3e1);
t281 = -0.8e1 * (0.3e1 / 0.4e1 * t275 + t255 + t379) * t170 - 0.8e1 * (t309 - t353) * t487 + 0.8e1 * (-pkin(2) * t184 + t183 * t252) * t278;
t280 = Icges(1,3) + Icges(2,3) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + Icges(3,1) / 0.2e1 + t458 + (0.2e1 * t244 ^ 2 + t256 + t259 + t378) * t464 + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2 + t278) * m(2);
t274 = pkin(3) * t275;
t254 = m(2) + m(3);
t210 = cos(t267);
t203 = sin(t269);
t202 = sin(t268);
t200 = sin(t266);
t199 = sin(t265);
t197 = sin(t263);
t196 = sin(t262);
t193 = t266 - pkin(7);
t192 = t266 + pkin(7);
t190 = -0.2e1 * pkin(7) + qJ(3,2);
t188 = -pkin(7) + t267;
t186 = pkin(7) + t267;
t185 = t261 + qJ(3,2);
t182 = t243 * pkin(1);
t180 = t241 * pkin(1);
t178 = t239 * pkin(1);
t165 = cos(t189);
t164 = cos(t188);
t163 = cos(t187);
t162 = cos(t186);
t154 = -Icges(3,4) + t364;
t147 = 0.2e1 * t189;
t146 = 0.2e1 * t187;
t132 = m(3) * t378 + Icges(3,3);
t125 = rSges(3,2) * t329;
t115 = (t272 / 0.2e1 - t271 / 0.2e1) * m(3) + t492;
t114 = g(1) * t176 - g(2) * t173;
t113 = g(1) * t175 - g(2) * t172;
t112 = g(1) * t174 - g(2) * t171;
t111 = g(1) * t173 + g(2) * t176;
t110 = g(1) * t172 + g(2) * t175;
t109 = g(1) * t171 + g(2) * t174;
t91 = t135 + t277 + 0.2e1 * t362 + t379;
t63 = -t176 * t350 + (pkin(3) * t403 - t499 * t176) * t242 + t144 * t403;
t62 = -t175 * t351 + (pkin(3) * t405 - t500 * t175) * t240 + t144 * t405;
t61 = -t174 * t352 + (pkin(3) * t406 - t501 * t174) * t238 + t144 * t406;
t60 = t173 * t350 + (pkin(3) * t399 + t499 * t173) * t242 + t144 * t399;
t59 = t172 * t351 + (pkin(3) * t401 + t500 * t172) * t240 + t144 * t401;
t58 = t171 * t352 + (pkin(3) * t402 + t501 * t171) * t238 + t144 * t402;
t56 = t213 * t456 + t128 * t327 + (-(-m(3) * t294 - t462) * t222 - t290) * t375 - t154 * t204 + t280;
t55 = t210 * t456 + t127 * t327 + (-(-m(3) * t295 - t462) * t222 - t290) * t375 - t154 * t201 + t280;
t54 = t207 * t456 + t126 * t327 + (-(-m(3) * t296 - t462) * t222 - t290) * t375 - t154 * t198 + t280;
t53 = -t121 * t254 * t334 + t348 * t452;
t52 = -t120 * t254 * t336 + t347 * t453;
t51 = -t119 * t254 * t337 + t349 * t454;
t50 = (-t176 * t317 + t254 * t63) * t414;
t49 = (-t175 * t318 + t254 * t62) * t416;
t48 = (-t174 * t319 + t254 * t61) * t417;
t47 = (t173 * t317 + t254 * t60) * t414;
t46 = (t172 * t318 + t254 * t59) * t416;
t45 = (t171 * t319 + t254 * t58) * t417;
t44 = t132 * t348 + (-t168 * t79 - t334 * t452) * t121;
t43 = t132 * t347 + (-t167 * t78 - t336 * t453) * t120;
t42 = t132 * t349 + (-t166 * t77 - t337 * t454) * t119;
t29 = (t79 * t420 + (t132 * t341 + t452 * t60) * t216) * t121;
t28 = (t78 * t421 + (t132 * t342 + t453 * t59) * t215) * t120;
t27 = (t77 * t427 + (t132 * t343 + t454 * t58) * t214) * t119;
t26 = (t79 * t428 + (-t132 * t338 + t452 * t63) * t216) * t121;
t25 = (t78 * t429 + (-t132 * t339 + t453 * t62) * t215) * t120;
t24 = (t77 * t430 + (-t132 * t340 + t454 * t61) * t214) * t119;
t21 = t477 * t67 - t359;
t20 = t477 * t66 - t360;
t19 = t477 * t65 - t361;
t18 = (t21 - t359) * t67 * t121;
t17 = (t20 - t360) * t66 * t120;
t16 = (t19 - t361) * t65 * t119;
t15 = ((-t144 - t181) * t468 * pkin(3) * t216 + (-(t219 * t275 + t91) * t67 + (-t144 * t242 * t67 + t434 * t477) * t373) * t67 * t396) * t121;
t14 = ((-t144 - t179) * t469 * pkin(3) * t215 + (-(t218 * t275 + t91) * t66 + (-t144 * t240 * t66 + t436 * t477) * t373) * t66 * t397) * t120;
t13 = ((-t144 - t177) * t470 * pkin(3) * t214 + (-(t217 * t275 + t91) * t65 + (-t144 * t238 * t65 + t437 * t477) * t373) * t65 * t398) * t119;
t12 = ((t212 * t358 + (t152 * t477 + t316 - t440) * t373) * t41 + (-t274 * t202 + t203 * t357 + t204 * t376 + t236 * t281) * t67) / (t418 + t385 / 0.2e1 - 0.2e1 * t354 + 0.2e1 * t288 + t308) * t67 * t455 + (t21 * t472 + t359 * t369 + (-0.2e1 * t440 + (-t212 + t242) * t323 + t475) * t67 + (-t275 * t202 + t203 * t371 + t204 * t467) * t41) / (0.4e1 * t288 + t310 - 0.4e1 * t354 + t385 + 0.2e1 * t418) * t41;
t11 = ((t206 * t358 + (t149 * t477 + t316 - t441) * t373) * t39 + (-t274 * t196 + t197 * t357 + t198 * t376 + t232 * t281) * t65) / (t419 + t387 / 0.2e1 - 0.2e1 * t355 + 0.2e1 * t289 + t308) * t65 * t455 + (t19 * t472 + t361 * t369 + (-0.2e1 * t441 + (-t206 + t238) * t323 + t475) * t65 + (-t275 * t196 + t197 * t371 + t198 * t467) * t39) / (0.4e1 * t289 + t310 - 0.4e1 * t355 + t387 + 0.2e1 * t419) * t39;
t10 = -m(2) * t111 - t254 * t15 + (-t128 * t12 - t303 * t468 - t111) * m(3);
t9 = -m(2) * t109 - t254 * t13 + (-t126 * t11 - t305 * t470 - t109) * m(3);
t8 = -t79 * t18 - t15 * t452 - t132 * t12 + 0.2e1 * (t154 * t219 + (t115 * t236 + t125) * t242 + t236 * t311 + t382) * t67 ^ 2 - m(3) * (t111 * t128 + (g(3) * t161 - t114 * t168) * t303);
t7 = -t77 * t16 - t13 * t454 - t132 * t11 + 0.2e1 * (t154 * t217 + (t115 * t232 + t125) * t238 + t232 * t311 + t382) * t65 ^ 2 - m(3) * (t109 * t126 + (g(3) * t159 - t112 * t166) * t305);
t6 = -t56 * t18 - t79 * t12 - 0.4e1 * ((t67 * t236 * t456 + t41 * t459) * t242 + t434 * t457 + (t219 - 0.1e1 / 0.2e1) * t67 * t154 + ((t242 * t297 - t306 * t434) * rSges(3,2) + (t242 * t306 * t41 + t236 * t297) * rSges(3,1)) * m(3)) * t41 - m(1) * (t114 * (-rSges(1,1) * t237 - rSges(1,2) * t243) - g(3) * (rSges(1,1) * t243 - rSges(1,2) * t237)) - m(2) * (t114 * (-rSges(2,1) * t161 - rSges(2,2) * t168 - t446) - g(3) * (rSges(2,1) * t168 - rSges(2,2) * t161 + t182)) - m(3) * (-t114 * t446 - g(3) * t182 + (-g(3) * t294 + t114 * t244) * t168 + (-t114 * t294 - t345) * t161);
t5 = -t54 * t16 - t77 * t11 - 0.4e1 * ((t65 * t232 * t456 + t39 * t459) * t238 + t437 * t457 + (t217 - 0.1e1 / 0.2e1) * t65 * t154 + ((t238 * t299 - t306 * t437) * rSges(3,2) + (t238 * t306 * t39 + t232 * t299) * rSges(3,1)) * m(3)) * t39 - m(1) * (t112 * (-rSges(1,1) * t233 - rSges(1,2) * t239) - g(3) * (rSges(1,1) * t239 - rSges(1,2) * t233)) - m(2) * (t112 * (-rSges(2,1) * t159 - rSges(2,2) * t166 - t448) - g(3) * (rSges(2,1) * t166 - rSges(2,2) * t159 + t178)) - m(3) * (-t112 * t448 - g(3) * t178 + (-g(3) * t296 + t112 * t244) * t166 + (-t112 * t296 - t345) * t159);
t4 = (t209 * t425 * t466 + t356 * t460 + (t450 + (-t170 - t179 / 0.2e1) * t252) * t377 + (-0.8e1 * (t275 / 0.4e1 + 0.3e1 / 0.2e1 * t277 + t381) * t449 - t274 * t199) * t66 - 0.4e1 * ((t283 + t493) * t158 + (t283 - t493) * t157) * pkin(1) - 0.3e1 * ((-t312 / 0.3e1 + (-t408 + t479) * t120 + t442) * sin(t193) + (t312 / 0.3e1 + (-t408 - t479) * t120 + t442) * sin(t192)) * pkin(1) * t275 + 0.4e1 * (cos(t190) - cos(t185)) * t252 * t64 + (sin(t146) * t370 + 0.4e1 * t164 * t363) * (t312 / 0.2e1 + (-t408 - t478) * t120 + t442) + (sin(t147) * t370 - 0.4e1 * t162 * t363) * (-t312 / 0.2e1 + (-t408 + t478) * t120 + t442) + (-0.12e2 * ((-t312 / 0.6e1 + (-t408 + t480) * t120 + t442) * sin(t188) + (t312 / 0.6e1 + (-t408 - t480) * t120 + t442) * sin(t186)) * t463 - 0.4e1 * (sin(t190) + sin(t185)) * t64 + 0.8e1 * (-t163 + t165) * pkin(1) * t426 + t66 * t200 * t465 + t252 * t210 * t377 - 0.16e2 * t309 * t435) * pkin(2)) / (t256 * t210 + t386 / 0.2e1 + (cos(t147) / 0.2e1 + cos(t146) / 0.2e1 + t210 - 0.1e1) * t278 + t380 * pkin(2) * t373 + ((t162 + t164 - 0.2e1 * t222) * t487 + (cos(t193) + cos(t192) - t163 - t165) * pkin(3)) * pkin(1) + t325) * t66 * t455 + (-t199 * t425 + (t360 + t426) * t369 + t20 * t472 + (t200 * t371 + t201 * t467) * t40 + (-0.2e1 * t89 * t210 - t380 * t323 + t356) * t66) / (0.2e1 * t102 * t210 + t386 + (t179 + t487) * t369 + 0.4e1 * (-t240 * pkin(2) + t144 * t209) * pkin(3) + t310) * t40;
t3 = -m(2) * t110 - t254 * t14 + (-t127 * t4 - t304 * t469 - t110) * m(3);
t2 = -t78 * t17 - t14 * t453 - t132 * t4 + 0.2e1 * (t154 * t218 + (t115 * t234 + t125) * t240 + t234 * t311 + t382) * t66 ^ 2 - m(3) * (t110 * t127 + (g(3) * t160 - t113 * t167) * t304);
t1 = -t55 * t17 - t78 * t4 - 0.4e1 * ((t40 * t459 + t435 * t456) * t240 + t436 * t457 + (t218 - 0.1e1 / 0.2e1) * t66 * t154 + ((t240 * t298 - t306 * t436) * rSges(3,2) + (t240 * t306 * t40 + t234 * t298) * rSges(3,1)) * m(3)) * t40 - m(1) * (t113 * (-rSges(1,1) * t235 - rSges(1,2) * t241) - g(3) * (rSges(1,1) * t241 - rSges(1,2) * t235)) - m(2) * (t113 * (-rSges(2,1) * t160 - rSges(2,2) * t167 - t447) - g(3) * (rSges(2,1) * t167 - rSges(2,2) * t160 + t180)) - m(3) * (-t113 * t447 - g(3) * t180 + (-g(3) * t295 + t113 * t244) * t167 + (-t113 * t295 - t345) * t160);
t22 = [(-g(1) + t231) * m(4) + (t24 * t445 + t25 * t444 + t26 * t443) * t391 + (t6 * t428 + ((-t26 * t338 + t50 * t63) * t231 + (t26 * t341 + t50 * t60) * t230 - t50 * t331 + t63 * t10 - t8 * t338) * t216 + t489 * (-t176 * t313 + t428 * t56) * t121) * t121 + (t1 * t429 + ((-t25 * t339 + t49 * t62) * t231 + (t25 * t342 + t49 * t59) * t230 - t49 * t332 + t62 * t3 - t2 * t339) * t215 + t490 * (-t175 * t314 + t429 * t55) * t120) * t120 + (t5 * t430 + ((-t24 * t340 + t48 * t61) * t231 + (t24 * t343 + t48 * t58) * t230 - t48 * t333 + t61 * t9 - t7 * t340) * t214 + (-t174 * t315 + t430 * t54) * t488) * t119; (-g(2) + t230) * m(4) + (t27 * t445 + t28 * t444 + t29 * t443) * t391 + (t6 * t420 + ((-t29 * t338 + t47 * t63) * t231 + (t29 * t341 + t47 * t60) * t230 - t47 * t331 + t60 * t10 + t8 * t341) * t216 + (t230 * t420 + t498) * (t173 * t313 + t420 * t56) * t121) * t121 + (t1 * t421 + ((-t28 * t339 + t46 * t62) * t231 + (t28 * t342 + t46 * t59) * t230 - t46 * t332 + t59 * t3 + t2 * t342) * t215 + (t230 * t421 + t497) * (t172 * t314 + t421 * t55) * t120) * t120 + (t5 * t427 + ((-t27 * t340 + t45 * t61) * t231 + (t27 * t343 + t45 * t58) * t230 - t45 * t333 + t58 * t9 + t7 * t343) * t214 + (t171 * t315 + t427 * t54) * t488) * t119; (-g(3) + t229) * m(4) + (t2 * t444 + t7 * t445 + t8 * t443 + (t42 * t445 + t43 * t444 + t44 * t443) * t229) * t276 + (-t168 * t6 + t489 * (-t121 * t168 * t56 + t348 * t79) + ((-t338 * t44 + t53 * t63) * t231 + (t341 * t44 + t53 * t60) * t230 + (-t229 * t53 - t10) * t431) * t216) * t121 + (-t167 * t1 + t490 * (-t120 * t167 * t55 + t347 * t78) + ((-t339 * t43 + t52 * t62) * t231 + (t342 * t43 + t52 * t59) * t230 + (-t229 * t52 - t3) * t432) * t215) * t120 + (-t166 * t5 + t491 * (-t119 * t166 * t54 + t349 * t77) + ((-t340 * t42 + t51 * t61) * t231 + (t343 * t42 + t51 * t58) * t230 + (-t229 * t51 - t9) * t433) * t214) * t119;];
tauX  = t22;
