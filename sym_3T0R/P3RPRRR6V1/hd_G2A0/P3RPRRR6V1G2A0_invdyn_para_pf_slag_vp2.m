% Calculate vector of inverse dynamics forces for parallel robot
% P3RPRRR6V1G2A0
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
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:37
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRRR6V1G2A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:36:06
% EndTime: 2020-08-06 18:36:17
% DurationCPUTime: 10.78s
% Computational Cost: add. (47637->651), mult. (53652->989), div. (8043->16), fcn. (36882->108), ass. (0->439)
t267 = xDDP(2);
t536 = t267 / 0.2e1;
t258 = cos(pkin(7));
t195 = t258 * pkin(1);
t157 = t195 + pkin(2);
t543 = mrSges(3,1) * t157;
t257 = sin(pkin(7));
t495 = pkin(1) * t257;
t350 = pkin(5) + t495;
t542 = t350 * mrSges(3,1);
t541 = t350 * mrSges(3,2);
t289 = -pkin(6) - pkin(5);
t421 = t257 * t289;
t138 = -pkin(1) + t421;
t270 = sin(qJ(1,3));
t118 = t138 * t270;
t276 = cos(qJ(1,3));
t416 = t276 * t289;
t424 = t257 * t276;
t540 = t118 - pkin(2) * t424 - (pkin(2) * t270 + t416) * t258;
t272 = sin(qJ(1,2));
t119 = t138 * t272;
t278 = cos(qJ(1,2));
t415 = t278 * t289;
t423 = t257 * t278;
t539 = t119 - pkin(2) * t423 - (pkin(2) * t272 + t415) * t258;
t274 = sin(qJ(1,1));
t120 = t138 * t274;
t280 = cos(qJ(1,1));
t414 = t280 * t289;
t422 = t257 * t280;
t538 = t120 - pkin(2) * t422 - (pkin(2) * t274 + t414) * t258;
t224 = qJ(1,3) + pkin(7);
t176 = sin(t224);
t266 = xDDP(3);
t263 = legFrame(3,2);
t145 = t263 + t224;
t146 = -t263 + t224;
t99 = -sin(t145) + sin(t146);
t537 = -t176 * t266 + t99 * t536;
t268 = xDDP(1);
t535 = t268 / 0.2e1;
t310 = pkin(3) ^ 2;
t311 = 0.1e1 / pkin(3);
t277 = cos(qJ(3,2));
t209 = t277 * pkin(3);
t166 = t209 + pkin(2);
t122 = 0.1e1 / (t195 + t166);
t271 = sin(qJ(3,2));
t251 = 0.1e1 / t271;
t446 = t122 * t251;
t443 = t166 * t272;
t74 = (t415 + t443) * t258 - t119 + t166 * t423;
t370 = t74 * t446;
t286 = xDP(3);
t217 = pkin(7) + qJ(3,2);
t190 = qJ(1,2) + t217;
t153 = cos(t190);
t221 = -pkin(7) + qJ(3,2);
t191 = qJ(1,2) - t221;
t154 = cos(t191);
t173 = sin(t217);
t175 = sin(t221);
t229 = qJ(1,2) + pkin(7);
t177 = sin(t229);
t186 = cos(t229);
t304 = 0.2e1 * qJ(3,2);
t236 = sin(t304);
t492 = pkin(3) * t236;
t493 = pkin(2) * t271;
t501 = 0.2e1 * t289;
t513 = -0.2e1 * pkin(2);
t514 = -0.2e1 * pkin(1);
t490 = (t177 * t501 + t278 * t514 + t186 * t513 + (-t153 - t154) * pkin(3)) / (0.2e1 * t493 + t492 + (t173 + t175) * pkin(1));
t371 = t286 * t490;
t264 = legFrame(2,2);
t205 = cos(t264);
t288 = xDP(1);
t429 = t205 * t288;
t202 = sin(t264);
t287 = xDP(2);
t433 = t202 * t287;
t533 = t310 * t311 * (t371 / 0.4e1 - (t429 / 0.4e1 - t433 / 0.4e1) * t370);
t147 = t264 + t229;
t148 = -t264 + t229;
t103 = cos(t148) + cos(t147);
t438 = t177 * t266;
t532 = t103 * t535 - t438;
t230 = qJ(1,1) + pkin(7);
t265 = legFrame(1,2);
t149 = t265 + t230;
t150 = -t265 + t230;
t104 = cos(t150) + cos(t149);
t178 = sin(t230);
t436 = t178 * t266;
t531 = t104 * t535 - t436;
t102 = cos(t146) + cos(t145);
t530 = t102 * t535 + t537;
t529 = 0.2e1 * pkin(2);
t528 = 0.4e1 * pkin(2);
t466 = t99 / 0.2e1;
t100 = -sin(t147) + sin(t148);
t456 = t100 / 0.2e1;
t101 = -sin(t149) + sin(t150);
t454 = t101 / 0.2e1;
t453 = t102 / 0.2e1;
t451 = t103 / 0.2e1;
t449 = t104 / 0.2e1;
t413 = t287 / 0.2e1;
t412 = t288 / 0.2e1;
t522 = t101 * t536;
t521 = t100 * t536;
t275 = cos(qJ(3,3));
t208 = t275 * pkin(3);
t165 = t208 + pkin(2);
t121 = 0.1e1 / (t195 + t165);
t201 = sin(t263);
t269 = sin(qJ(3,3));
t250 = 0.1e1 / t269;
t444 = t165 * t270;
t463 = t311 * ((t416 + t444) * t258 - t118 + t165 * t424);
t362 = t250 * t463;
t337 = t201 * t362;
t215 = pkin(7) + qJ(3,3);
t188 = qJ(1,3) + t215;
t151 = cos(t188);
t219 = -pkin(7) + qJ(3,3);
t189 = qJ(1,3) - t219;
t152 = cos(t189);
t172 = sin(t215);
t174 = sin(t219);
t185 = cos(t224);
t301 = 0.2e1 * qJ(3,3);
t233 = sin(t301);
t491 = (t176 * t501 + t276 * t514 + t185 * t513 + (-t151 - t152) * pkin(3)) / (t269 * t529 + pkin(3) * t233 + (t172 + t174) * pkin(1));
t374 = t311 * t491;
t487 = t121 * t287 * t337 + t286 * t374;
t462 = t311 * t74;
t360 = t251 * t462;
t520 = (t429 / 0.6e1 - t433 / 0.6e1) * t360;
t519 = (t429 / 0.3e1 - t433 / 0.3e1) * t360;
t518 = (t429 / 0.2e1 - t433 / 0.2e1) * t360;
t298 = 0.2e1 * pkin(7);
t213 = cos(t298);
t313 = pkin(1) ^ 2;
t171 = t313 * t213;
t517 = -t171 - t310 / 0.2e1;
t511 = -0.4e1 * pkin(3);
t510 = 2 * g(3);
t204 = cos(t263);
t336 = t204 * t362;
t331 = t288 * t336;
t32 = -t121 * t331 + t487;
t509 = t32 ^ 2;
t33 = (t371 + (-t429 + t433) * t370) * t311;
t508 = t33 ^ 2;
t203 = sin(t265);
t206 = cos(t265);
t279 = cos(qJ(3,1));
t210 = t279 * pkin(3);
t168 = t210 + pkin(2);
t123 = 0.1e1 / (t195 + t168);
t273 = sin(qJ(3,1));
t252 = 0.1e1 / t273;
t445 = t123 * t252;
t388 = pkin(7) + qJ(3,1);
t192 = qJ(1,1) + t388;
t155 = cos(t192);
t389 = -pkin(7) + qJ(3,1);
t193 = qJ(1,1) - t389;
t156 = cos(t193);
t187 = cos(t230);
t307 = 0.2e1 * qJ(3,1);
t239 = sin(t307);
t489 = (t178 * t501 + t280 * t514 + t187 * t513 + (-t155 - t156) * pkin(3)) / (t273 * t529 + pkin(3) * t239 + (sin(t388) + sin(t389)) * pkin(1));
t442 = t168 * t274;
t76 = (t414 + t442) * t258 - t120 + t168 * t422;
t34 = (t286 * t489 + (t203 * t287 - t206 * t288) * t76 * t445) * t311;
t507 = t34 ^ 2;
t506 = -2 * Ifges(3,4);
t386 = pkin(2) * t195;
t312 = pkin(2) ^ 2;
t293 = 0.2e1 * t312;
t402 = t293 + t313;
t105 = 0.4e1 * t386 + t171 + t402;
t505 = -0.2e1 * t310 - 0.2e1 * t105;
t170 = m(3) * pkin(5) - mrSges(2,2) + mrSges(3,3);
t504 = -0.2e1 * t170;
t207 = m(3) * pkin(2) + mrSges(2,1);
t503 = -0.2e1 * t207;
t502 = -0.2e1 * t258;
t500 = -0.6e1 * t310;
t297 = -0.2e1 * t313;
t499 = pkin(2) * pkin(3);
t291 = m(2) + m(3);
t498 = mrSges(3,2) * pkin(2);
t497 = pkin(3) * t33;
t496 = t311 / 0.2e1;
t143 = pkin(6) + t350;
t494 = pkin(2) * t143;
t488 = 0.2e1 * t487;
t486 = (t102 * t412 + t413 * t99) * t121;
t85 = t100 * t122 * t413;
t87 = t103 * t122 * t412;
t485 = t85 + t87;
t115 = g(1) * t204 - g(2) * t201;
t484 = mrSges(3,1) * t115;
t116 = g(1) * t205 - g(2) * t202;
t483 = mrSges(3,1) * t116;
t117 = g(1) * t206 - g(2) * t203;
t482 = mrSges(3,1) * t117;
t481 = mrSges(3,1) * t269;
t480 = mrSges(3,1) * t271;
t479 = mrSges(3,1) * t273;
t478 = mrSges(3,2) * t115;
t477 = mrSges(3,2) * t116;
t476 = mrSges(3,2) * t117;
t475 = mrSges(3,2) * t273;
t248 = cos(t307);
t390 = -0.2e1 * t195;
t212 = sin(t298);
t407 = t313 * t212;
t341 = t289 * t390 + t407;
t95 = t341 + 0.2e1 * t494;
t474 = t248 * t95;
t472 = t269 * t32;
t471 = t271 * t33;
t470 = t273 * t34;
t73 = (t165 * t276 - t270 * t289) * t258 - t138 * t276 - t257 * t444;
t469 = t275 * t73;
t75 = (t166 * t278 - t272 * t289) * t258 - t138 * t278 - t257 * t443;
t468 = t277 * t75;
t77 = (t168 * t280 - t274 * t289) * t258 - t138 * t280 - t257 * t442;
t467 = t279 * t77;
t437 = t177 * t286;
t353 = t122 * t437;
t65 = -t353 + t485;
t465 = t289 * t65;
t464 = t310 * t33;
t461 = t311 * t76;
t63 = t313 * t65;
t253 = t275 ^ 2;
t395 = 0.2e1 * pkin(3);
t427 = t250 * t275;
t439 = t176 * t286;
t64 = -t121 * t439 + t486;
t378 = pkin(1) * t421;
t140 = -0.2e1 * t378;
t256 = t289 ^ 2;
t403 = t256 + t313;
t97 = t140 + t312 + 0.2e1 * t386 + t403;
t460 = -((-t157 - t208) * t509 * pkin(3) * t250 + (-(t253 * t310 + t97) * t64 + (-t157 * t275 * t64 + t143 * t472) * t395) * t64 * t427) * t121 - g(1) * t201 - g(2) * t204;
t254 = t277 ^ 2;
t426 = t251 * t277;
t459 = -((-t157 - t209) * t508 * pkin(3) * t251 + (-(t254 * t310 + t97) * t65 + (-t157 * t277 * t65 + t143 * t471) * t395) * t65 * t426) * t122 - g(1) * t202 - g(2) * t205;
t255 = t279 ^ 2;
t425 = t252 * t279;
t66 = (t101 * t413 + t104 * t412 - t178 * t286) * t123;
t458 = ((-t157 - t210) * t507 * pkin(3) * t252 + (-(t255 * t310 + t97) * t66 + (-t157 * t279 * t66 + t143 * t470) * t395) * t66 * t425) * t123 + g(1) * t203 + g(2) * t206;
t448 = t105 * t248;
t447 = t121 * t250;
t441 = t170 * t257;
t435 = t201 * t269;
t434 = t202 * t271;
t432 = t203 * t273;
t431 = t204 * t269;
t430 = t205 * t271;
t428 = t206 * t273;
t262 = Ifges(3,1) - Ifges(3,2);
t420 = t262 * t269;
t419 = t262 * t271;
t418 = t262 * t273;
t417 = t266 * t311;
t299 = 0.4e1 * qJ(3,3);
t411 = t310 * cos(t299);
t302 = 0.4e1 * qJ(3,2);
t410 = t310 * cos(t302);
t305 = 0.4e1 * qJ(3,1);
t409 = t310 * cos(t305);
t405 = t256 / 0.2e1 + t313;
t303 = 0.3e1 * qJ(3,2);
t244 = cos(t303);
t404 = t244 - t277;
t401 = -0.2e1 * t498;
t292 = 0.3e1 * t312;
t296 = 0.2e1 * t313;
t399 = (0.6e1 * t386 + t140 + t296 + t292 + t256 - t517) * t511;
t398 = 0.2e1 * pkin(1);
t393 = t157 * t511;
t392 = pkin(3) * t297;
t391 = -0.4e1 * t195;
t387 = pkin(1) * pkin(3) * t289;
t385 = pkin(3) * t472;
t384 = pkin(3) * t471;
t383 = pkin(3) * t470;
t382 = -0.2e1 * t143 * t310;
t381 = t157 * t500;
t380 = 0.2e1 * t407;
t167 = t210 + t529;
t379 = t167 * t195;
t377 = pkin(3) * (t258 * t270 + t424) * t253;
t376 = pkin(3) * (t258 * t272 + t423) * t254;
t375 = pkin(3) * (t258 * t274 + t422) * t255;
t373 = t311 * t489;
t372 = t311 * t490;
t369 = t201 * t463;
t368 = t202 * t462;
t367 = t203 * t461;
t366 = t204 * t463;
t365 = t205 * t462;
t364 = t206 * t461;
t363 = t73 * t427;
t361 = t75 * t426;
t359 = t77 * t425;
t358 = t266 * t469;
t357 = t266 * t468;
t356 = t266 * t467;
t214 = pkin(7) + t301;
t218 = -pkin(7) + t301;
t352 = cos(t214) + cos(t218) + t502;
t351 = -0.2e1 * t312 + t517;
t300 = 0.3e1 * qJ(3,3);
t241 = cos(t300);
t348 = pkin(3) * (t241 - t275);
t347 = t143 * t395;
t345 = t253 * t506 + Ifges(3,4);
t344 = t254 * t506 + Ifges(3,4);
t343 = t255 * t506 + Ifges(3,4);
t342 = -t289 + t495;
t128 = Ifges(3,6) - t541;
t326 = Ifges(3,5) - t542;
t89 = t128 * t277 + t271 * t326;
t340 = t89 * t360;
t90 = t128 * t279 + t273 * t326;
t339 = t252 * t90 * t461;
t338 = t311 * t371;
t335 = -0.2e1 * t171 - 0.4e1 * t312 + t297 - t310;
t334 = 0.3e1 / 0.8e1 * t310 + t312 / 0.2e1 + t405;
t333 = -t313 + t351;
t329 = -Ifges(3,5) / 0.2e1 + t542 / 0.2e1;
t328 = Ifges(3,6) / 0.2e1 - t541 / 0.2e1;
t327 = mrSges(3,2) * t195 + t498;
t306 = 0.3e1 * qJ(3,1);
t247 = cos(t306);
t322 = pkin(3) * (-pkin(2) * t279 + t157 * t247);
t317 = t63 + (0.3e1 / 0.4e1 * t87 + 0.3e1 / 0.4e1 * t85 - 0.3e1 / 0.4e1 * t353) * t310 + 0.3e1 * (t312 + t256 / 0.3e1) * t65;
t316 = -0.8e1 * (0.3e1 / 0.4e1 * t310 + t292 + t403) * t195 - 0.8e1 * (t334 - t378) * t529 + 0.8e1 * (-pkin(2) * t213 + t212 * t289) * t313;
t315 = (pkin(5) ^ 2 + t312) * m(3) + 0.2e1 * mrSges(3,3) * pkin(5) + t291 * t313 + Ifges(3,1) + Ifges(1,3) + Ifges(2,3);
t309 = pkin(3) * t310;
t284 = mrSges(3,1) * g(3);
t283 = mrSges(1,2) * g(3);
t282 = mrSges(3,2) * g(3);
t245 = cos(t304);
t242 = cos(t301);
t238 = sin(t306);
t237 = sin(t305);
t235 = sin(t303);
t234 = sin(t302);
t232 = sin(t300);
t231 = sin(t299);
t228 = qJ(3,2) - 0.2e1 * pkin(7);
t227 = qJ(3,2) + t298;
t226 = t303 - pkin(7);
t225 = t303 + pkin(7);
t223 = t300 - pkin(7);
t222 = t300 + pkin(7);
t220 = -pkin(7) + t304;
t216 = pkin(7) + t304;
t199 = t271 * mrSges(3,2);
t198 = t269 * mrSges(3,2);
t184 = cos(t221);
t183 = cos(t220);
t181 = cos(t217);
t180 = cos(t216);
t164 = t208 + t529;
t162 = 0.2e1 * t221;
t161 = 0.2e1 * t219;
t160 = 0.2e1 * t217;
t159 = 0.2e1 * t215;
t144 = t291 * pkin(1) + mrSges(1,1);
t142 = t207 * t510;
t137 = mrSges(3,2) * t157;
t136 = t144 * g(3);
t135 = t170 * t510;
t134 = mrSges(3,2) * t279 + t479;
t133 = mrSges(3,2) * t277 + t480;
t132 = mrSges(3,2) * t275 + t481;
t131 = mrSges(3,1) * t277 - t199;
t130 = mrSges(3,1) * t275 - t198;
t129 = -t279 * mrSges(3,1) + t475;
t88 = t128 * t275 + t269 * t326;
t69 = -t262 * t255 + 0.2e1 * (Ifges(3,4) * t273 + t543) * t279 + t273 * t401 + (-(-t207 + t475) * t258 + t441) * t398 + t315;
t68 = -t262 * t254 + 0.2e1 * (Ifges(3,4) * t271 + t543) * t277 + t271 * t401 + (-(t199 - t207) * t258 + t441) * t398 + t315;
t67 = -t262 * t253 + 0.2e1 * (Ifges(3,4) * t269 + t543) * t275 + t269 * t401 + (-(t198 - t207) * t258 + t441) * t398 + t315;
t62 = t206 * t375 + (pkin(3) * t432 - t538 * t206) * t279 + t157 * t432;
t61 = -t203 * t375 + (pkin(3) * t428 + t538 * t203) * t279 + t157 * t428;
t60 = t205 * t376 + (pkin(3) * t434 - t539 * t205) * t277 + t157 * t434;
t59 = -t202 * t376 + (pkin(3) * t430 + t539 * t202) * t277 + t157 * t430;
t58 = t204 * t377 + (pkin(3) * t435 - t540 * t204) * t275 + t157 * t435;
t57 = -t201 * t377 + (pkin(3) * t431 + t540 * t201) * t275 + t157 * t431;
t55 = t123 * t291 * t359 - t129 * t373;
t54 = t122 * t291 * t361 + t131 * t372;
t53 = t121 * t291 * t363 + t130 * t374;
t49 = (t129 * t364 + t291 * t62) * t445;
t48 = (-t129 * t367 + t291 * t61) * t445;
t47 = (-t131 * t365 + t291 * t60) * t446;
t46 = (t131 * t368 + t291 * t59) * t446;
t45 = (-t130 * t366 + t291 * t58) * t447;
t44 = (t130 * t369 + t291 * t57) * t447;
t43 = Ifges(3,3) * t373 + (-t129 * t359 - t178 * t90) * t123;
t42 = Ifges(3,3) * t372 + (t131 * t361 - t177 * t89) * t122;
t41 = Ifges(3,3) * t374 + (t130 * t363 - t176 * t88) * t121;
t31 = (t90 * t449 + (-Ifges(3,3) * t364 - t129 * t62) * t252) * t123;
t30 = (t89 * t451 + (-Ifges(3,3) * t365 + t131 * t60) * t251) * t122;
t29 = (t88 * t453 + (-Ifges(3,3) * t366 + t130 * t58) * t250) * t121;
t28 = (t90 * t454 + (Ifges(3,3) * t367 - t129 * t61) * t252) * t123;
t27 = (t89 * t456 + (Ifges(3,3) * t368 + t131 * t59) * t251) * t122;
t26 = (t88 * t466 + (Ifges(3,3) * t369 + t130 * t57) * t250) * t121;
t25 = (-0.2e1 * t331 - t439) * t121 + t486 + t488;
t24 = (0.2e1 * t331 - t439) * t121 + t486 - t488;
t21 = t342 * t65 - t384;
t20 = t342 * t64 - t385;
t19 = t342 * t66 - t383;
t18 = (t21 - t384) * t65 * t122;
t17 = (t20 - t385) * t64 * t121;
t16 = (t19 - t383) * t66 * t123;
t12 = ((t247 * t382 + (t143 * t167 + t341 - t474) * t395) * t34 + (-t309 * t237 + t238 * t381 + t239 * t399 + t316 * t273) * t66) / (t448 + t409 / 0.2e1 - 0.2e1 * t379 + 0.2e1 * t322 + t333) * t66 * t496 + (t19 * t528 + t383 * t391 + (-0.2e1 * t474 + t380 + t289 * t391 + (-t247 + t279) * t347) * t66 + (-t310 * t237 + t238 * t393 + t239 * t505) * t34) / (0.4e1 * t322 + t335 - 0.4e1 * t379 + t409 + 0.2e1 * t448) * t34;
t11 = t12 * t129 - t134 * t507 - t458 * t291;
t10 = ((t241 * t382 + (t143 * t164 - t242 * t95 + t341) * t395) * t32 + (-t309 * t231 + t232 * t381 + t233 * t399 + t316 * t269) * t64) / (t105 * t242 + t411 / 0.2e1 + t164 * t390 + (-t275 * pkin(2) + t157 * t241) * t395 + t333) * t64 * t496 + (t20 * t528 + (((t331 - t439) * t121 + t486 - t487) * sin(t161) - ((-t331 - t439) * t121 + t486 + t487) * sin(t159)) * t313 + (-0.2e1 * (t310 + t402) * t233 - 0.4e1 * t232 * t499 - t310 * t231) * t32 + (t380 + (t242 * t528 + 0.2e1 * t348) * t289) * t64 + ((sin(t218) * t24 - sin(t214) * t25) * t529 + t352 * t64 * t501 + ((-sin(t222) - t174) * t25 + (sin(t223) + t172) * t24) * pkin(3)) * pkin(1)) / ((t296 + 0.4e1 * t312) * t242 + t411 + (cos(t161) + cos(t159)) * t313 + t348 * t528 + (t352 * t528 + (cos(t223) + cos(t222) - cos(t219) - cos(t215)) * t395) * pkin(1) + t335) * t32;
t9 = -t90 * t16 - Ifges(3,3) * t12 + t66 ^ 2 * ((t137 - t418) * t279 + t157 * t479 + t343) + (g(3) * t187 + t117 * t178) * t134 + t458 * t129;
t8 = -t69 * t16 - t90 * t12 - 0.2e1 * ((t328 * t273 + t329 * t279) * t34 + ((t327 - t418) * t279 + t273 * t543 + t343) * t66) * t34 + (-t282 - t482) * t156 / 0.2e1 + (t284 - t476) * sin(t193) / 0.2e1 + (t282 - t482) * t155 / 0.2e1 + (t284 + t476) * sin(t192) / 0.2e1 + (t117 * t503 - t135) * t187 / 0.2e1 + (t117 * t504 + t142) * t178 / 0.2e1 + (-t117 * t144 + t283) * t280 + (mrSges(1,2) * t117 + t136) * t274;
t7 = (t244 * t464 * t501 + t380 * t497 + 0.4e1 * (t494 + (pkin(2) * t245 - t195 - t209 / 0.2e1) * t289) * t497 + (sin(t160) * t392 + 0.4e1 * t183 * t387) * (t338 / 0.2e1 + (-t437 - t518) * t122 + t485) + (sin(t162) * t392 - 0.4e1 * t180 * t387) * (-t338 / 0.2e1 + (-t437 + t518) * t122 + t485) + (-0.8e1 * (t310 / 0.4e1 + 0.3e1 / 0.2e1 * t312 + t405) * t492 + pkin(2) * t235 * t500 - t309 * t234 - 0.16e2 * t334 * t493) * t65 + 0.4e1 * ((cos(t228) - cos(t227)) * t289 - (sin(t228) + sin(t227)) * pkin(2)) * t63 + (-0.4e1 * (t317 + t533) * t175 - 0.4e1 * (t317 - t533) * t173 - 0.3e1 * ((-t338 / 0.3e1 + (-t437 + t519) * t122 + t485) * sin(t226) + (t338 / 0.3e1 + (-t437 - t519) * t122 + t485) * sin(t225)) * t310 - 0.12e2 * ((-t338 / 0.6e1 + (-t437 + t520) * t122 + t485) * sin(t220) + (t338 / 0.6e1 + (-t437 - t520) * t122 + t485) * sin(t216)) * t499 + 0.8e1 * (-t181 + t184) * pkin(2) * t465) * pkin(1)) / (t293 * t245 + t410 / 0.2e1 + (cos(t162) / 0.2e1 + cos(t160) / 0.2e1 + t245 - 0.1e1) * t313 + t404 * pkin(2) * t395 + ((t180 + t183 + t502) * t529 + (cos(t226) + cos(t225) - t181 - t184) * pkin(3)) * pkin(1) + t351) * t65 * t496 + (-t234 * t464 + (t384 + t465) * t391 + t21 * t528 + (t235 * t393 + t236 * t505) * t33 + (-0.2e1 * t95 * t245 - t404 * t347 + t380) * t65) / (0.2e1 * t105 * t245 + t410 + (t209 + t529) * t391 + 0.4e1 * (-t277 * pkin(2) + t157 * t244) * pkin(3) + t335) * t33;
t6 = -t10 * t130 - t132 * t509 + t460 * t291;
t5 = -t88 * t17 - Ifges(3,3) * t10 + t64 ^ 2 * ((t137 - t420) * t275 + t157 * t481 + t345) + (g(3) * t185 + t115 * t176) * t132 + t460 * t130;
t4 = -t67 * t17 - t88 * t10 - 0.2e1 * ((t328 * t269 + t329 * t275) * t32 + ((t327 - t420) * t275 + t269 * t543 + t345) * t64) * t32 + (-t282 - t484) * t152 / 0.2e1 + (t284 - t478) * sin(t189) / 0.2e1 + (t282 - t484) * t151 / 0.2e1 + (t284 + t478) * sin(t188) / 0.2e1 + (t115 * t503 - t135) * t185 / 0.2e1 + (t115 * t504 + t142) * t176 / 0.2e1 + (-t115 * t144 + t283) * t276 + (mrSges(1,2) * t115 + t136) * t270;
t3 = -t131 * t7 - t133 * t508 + t459 * t291;
t2 = -t89 * t18 - Ifges(3,3) * t7 + t65 ^ 2 * ((t137 - t419) * t277 + t157 * t480 + t344) + (g(3) * t186 + t116 * t177) * t133 + t459 * t131;
t1 = -t68 * t18 - t89 * t7 - 0.2e1 * ((t328 * t271 + t329 * t277) * t33 + ((t327 - t419) * t277 + t271 * t543 + t344) * t65) * t33 + (-t282 - t483) * t154 / 0.2e1 + (t284 - t477) * sin(t191) / 0.2e1 + (t282 - t483) * t153 / 0.2e1 + (t284 + t477) * sin(t190) / 0.2e1 + (t116 * t503 - t135) * t186 / 0.2e1 + (t116 * t504 + t142) * t177 / 0.2e1 + (-t116 * t144 + t283) * t278 + (mrSges(1,2) * t116 + t136) * t272;
t13 = [(-g(1) + t268) * m(4) + (t29 * t491 + t30 * t490 + t31 * t489) * t417 + (t8 * t449 + ((-t31 * t364 + t49 * t62) * t268 + (t31 * t367 + t49 * t61) * t267 + t49 * t356 + t62 * t11 - t9 * t364) * t252 + (t268 * t449 - t436 + t522) * (-t206 * t339 + t69 * t449) * t123) * t123 + (t1 * t451 + ((-t30 * t365 + t47 * t60) * t268 + (t30 * t368 + t47 * t59) * t267 + t47 * t357 + t60 * t3 - t2 * t365) * t251 + (t268 * t451 - t438 + t521) * (-t205 * t340 + t68 * t451) * t122) * t122 + (t4 * t453 + ((-t29 * t366 + t45 * t58) * t268 + (t29 * t369 + t45 * t57) * t267 + t45 * t358 + t58 * t6 - t5 * t366) * t250 + (t268 * t453 + t537) * (-t88 * t336 + t67 * t453) * t121) * t121; (-g(2) + t267) * m(4) + (t26 * t491 + t27 * t490 + t28 * t489) * t417 + (t8 * t454 + ((-t28 * t364 + t48 * t62) * t268 + (t28 * t367 + t48 * t61) * t267 + t48 * t356 + t61 * t11 + t9 * t367) * t252 + (t267 * t454 + t531) * (t203 * t339 + t69 * t454) * t123) * t123 + (t1 * t456 + ((-t27 * t365 + t46 * t60) * t268 + (t27 * t368 + t46 * t59) * t267 + t46 * t357 + t59 * t3 + t2 * t368) * t251 + (t267 * t456 + t532) * (t202 * t340 + t68 * t456) * t122) * t122 + (t4 * t466 + ((-t26 * t366 + t44 * t58) * t268 + (t26 * t369 + t44 * t57) * t267 + t44 * t358 + t57 * t6 + t5 * t369) * t250 + t530 * (t88 * t337 + t67 * t466) * t121) * t121; (-g(3) + t266) * m(4) + (t2 * t490 + t5 * t491 + t9 * t489 + (t41 * t491 + t42 * t490 + t43 * t489) * t266) * t311 + (-t178 * t8 + (t522 + t531) * (-t123 * t178 * t69 + t90 * t373) + ((-t43 * t364 + t55 * t62) * t268 + (t43 * t367 + t55 * t61) * t267 + (t266 * t55 + t11) * t467) * t252) * t123 + (-t177 * t1 + (t521 + t532) * (-t122 * t177 * t68 + t89 * t372) + ((-t42 * t365 + t54 * t60) * t268 + (t42 * t368 + t54 * t59) * t267 + (t266 * t54 + t3) * t468) * t251) * t122 + (-t176 * t4 + t530 * (-t121 * t176 * t67 + t88 * t374) + ((-t41 * t366 + t53 * t58) * t268 + (t41 * t369 + t53 * t57) * t267 + (t266 * t53 + t6) * t469) * t250) * t121;];
tauX  = t13;
