% Calculate vector of inverse dynamics forces for parallel robot
% P3RRPRR8V2G1A0
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2020-08-06 21:04
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:02:41
% EndTime: 2020-08-06 21:02:51
% DurationCPUTime: 9.70s
% Computational Cost: add. (31146->587), mult. (38811->920), div. (5487->15), fcn. (28305->62), ass. (0->393)
t235 = -qJ(3,1) - pkin(5);
t144 = mrSges(3,1) * t235 + Ifges(3,5);
t227 = sin(pkin(7));
t228 = cos(pkin(7));
t458 = mrSges(3,2) * t235 + Ifges(3,6);
t483 = -t144 * t228 + t458 * t227;
t234 = -qJ(3,2) - pkin(5);
t143 = mrSges(3,1) * t234 + Ifges(3,5);
t459 = mrSges(3,2) * t234 + Ifges(3,6);
t482 = -t143 * t228 + t459 * t227;
t233 = -qJ(3,3) - pkin(5);
t142 = mrSges(3,1) * t233 + Ifges(3,5);
t460 = mrSges(3,2) * t233 + Ifges(3,6);
t481 = -t142 * t228 + t460 * t227;
t480 = (m(2) * pkin(5));
t427 = pkin(3) * t228;
t152 = pkin(2) + t427;
t251 = cos(qJ(2,3));
t479 = t152 * t251;
t253 = cos(qJ(2,2));
t478 = t152 * t253;
t255 = cos(qJ(2,1));
t477 = t152 * t255;
t205 = -pkin(6) + t235;
t194 = 0.1e1 / t205;
t258 = xDP(3);
t201 = t255 * pkin(2);
t215 = qJ(2,1) + pkin(7);
t179 = cos(t215);
t429 = pkin(3) * t179;
t127 = t201 + t429;
t118 = 0.1e1 / t127;
t249 = sin(qJ(2,1));
t124 = pkin(2) * t249 + pkin(3) * sin(t215);
t275 = 0.2e1 * qJ(2,1);
t214 = t275 + pkin(7);
t167 = sin(t214);
t475 = 0.2e1 * pkin(2);
t366 = pkin(3) * t475;
t280 = pkin(2) ^ 2;
t368 = t280 * sin(t275);
t155 = 0.2e1 * t215;
t278 = pkin(3) ^ 2;
t371 = t278 * sin(t155);
t293 = t167 * t366 + t368 + t371;
t476 = 2 * pkin(1);
t75 = t124 * t476 + t293;
t442 = t75 / 0.2e1;
t334 = t118 * t442;
t260 = xDP(1);
t232 = legFrame(1,3);
t185 = sin(t232);
t188 = cos(t232);
t158 = t201 + pkin(1);
t250 = sin(qJ(1,1));
t256 = cos(qJ(1,1));
t89 = t158 * t250 + t205 * t256;
t92 = t158 * t256 - t205 * t250;
t96 = -t185 * t250 + t188 * t256;
t63 = -t185 * t89 + t188 * t92 + t96 * t429;
t401 = t260 * t63;
t259 = xDP(2);
t99 = t185 * t256 + t188 * t250;
t66 = t185 * t92 + t188 * t89 + t99 * t429;
t404 = t259 * t66;
t33 = (t258 * t334 + t401 + t404) * t194;
t204 = -pkin(6) + t234;
t193 = 0.1e1 / t204;
t200 = t253 * pkin(2);
t213 = qJ(2,2) + pkin(7);
t175 = cos(t213);
t430 = pkin(3) * t175;
t126 = t200 + t430;
t115 = 0.1e1 / t126;
t247 = sin(qJ(2,2));
t123 = pkin(2) * t247 + pkin(3) * sin(t213);
t273 = 0.2e1 * qJ(2,2);
t212 = t273 + pkin(7);
t165 = sin(t212);
t369 = t280 * sin(t273);
t154 = 0.2e1 * t213;
t372 = t278 * sin(t154);
t294 = t165 * t366 + t369 + t372;
t74 = t123 * t476 + t294;
t443 = t74 / 0.2e1;
t335 = t115 * t443;
t231 = legFrame(2,3);
t184 = sin(t231);
t187 = cos(t231);
t157 = t200 + pkin(1);
t248 = sin(qJ(1,2));
t254 = cos(qJ(1,2));
t88 = t157 * t248 + t204 * t254;
t91 = t157 * t254 - t204 * t248;
t95 = -t184 * t248 + t187 * t254;
t62 = -t184 * t88 + t187 * t91 + t95 * t430;
t402 = t260 * t62;
t98 = t184 * t254 + t187 * t248;
t65 = t184 * t91 + t187 * t88 + t98 * t430;
t405 = t259 * t65;
t32 = (t258 * t335 + t402 + t405) * t193;
t203 = -pkin(6) + t233;
t192 = 0.1e1 / t203;
t199 = t251 * pkin(2);
t211 = qJ(2,3) + pkin(7);
t172 = cos(t211);
t431 = pkin(3) * t172;
t125 = t199 + t431;
t112 = 0.1e1 / t125;
t245 = sin(qJ(2,3));
t122 = pkin(2) * t245 + pkin(3) * sin(t211);
t271 = 0.2e1 * qJ(2,3);
t210 = t271 + pkin(7);
t163 = sin(t210);
t370 = t280 * sin(t271);
t153 = 0.2e1 * t211;
t373 = t278 * sin(t153);
t295 = t163 * t366 + t370 + t373;
t73 = t122 * t476 + t295;
t444 = t73 / 0.2e1;
t336 = t112 * t444;
t230 = legFrame(3,3);
t183 = sin(t230);
t186 = cos(t230);
t156 = t199 + pkin(1);
t246 = sin(qJ(1,3));
t252 = cos(qJ(1,3));
t87 = t156 * t246 + t203 * t252;
t90 = t156 * t252 - t203 * t246;
t94 = -t183 * t246 + t186 * t252;
t61 = -t183 * t87 + t186 * t90 + t94 * t431;
t403 = t260 * t61;
t97 = t183 * t252 + t186 * t246;
t64 = t183 * t90 + t186 * t87 + t97 * t431;
t406 = t259 * t64;
t31 = (t258 * t336 + t403 + t406) * t192;
t474 = -t258 / 0.2e1;
t113 = 0.1e1 / t125 ^ 2;
t226 = t258 ^ 2;
t473 = t113 * t226;
t116 = 0.1e1 / t126 ^ 2;
t472 = t116 * t226;
t119 = 0.1e1 / t127 ^ 2;
t471 = t119 * t226;
t209 = t228 ^ 2;
t467 = (-t209 + 0.1e1) * t278;
t181 = t227 * mrSges(3,1);
t239 = Ifges(3,2) - Ifges(3,1);
t80 = (-pkin(2) * mrSges(3,2) - t227 * t239) * t228 - pkin(2) * t181 + 0.2e1 * Ifges(3,4) * t209 + Ifges(2,4) - Ifges(3,4);
t189 = (t233 * m(3));
t241 = (mrSges(2,3) + mrSges(3,3));
t324 = -mrSges(1,2) + t241 + t480;
t466 = -t324 + t189;
t190 = (t234 * m(3));
t465 = -t324 + t190;
t191 = (t235 * m(3));
t464 = -t324 + t191;
t375 = t227 * t249;
t346 = pkin(3) * t375;
t463 = -t346 + t477;
t376 = t227 * t247;
t347 = pkin(3) * t376;
t462 = -t347 + t478;
t377 = t227 * t245;
t348 = pkin(3) * t377;
t461 = -t348 + t479;
t456 = 0.2e1 * mrSges(3,1);
t455 = 2 * mrSges(3,3);
t81 = 0.1e1 / t461;
t428 = pkin(3) * t227;
t84 = t152 * t245 + t251 * t428;
t419 = t84 * t81;
t55 = (t258 * t419 + t259 * t97 + t260 * t94) * t192;
t454 = -0.2e1 * t55;
t82 = 0.1e1 / t462;
t85 = t152 * t247 + t253 * t428;
t418 = t85 * t82;
t56 = (t258 * t418 + t259 * t98 + t260 * t95) * t193;
t453 = -0.2e1 * t56;
t83 = 0.1e1 / t463;
t86 = t152 * t249 + t255 * t428;
t417 = t86 * t83;
t57 = (t258 * t417 + t259 * t99 + t260 * t96) * t194;
t452 = -0.2e1 * t57;
t162 = pkin(2) * t427;
t365 = t162 + t280 / 0.2e1;
t451 = -0.2e1 * (t209 - 0.1e1 / 0.2e1) * t278 - 0.2e1 * t365;
t450 = -0.4e1 * pkin(1) * (t278 / 0.2e1 + t365);
t279 = pkin(2) * t280;
t432 = pkin(2) * t278;
t449 = -0.2e1 * t279 - 0.4e1 * t432;
t448 = 0.2e1 * t228;
t265 = m(3) * pkin(2);
t447 = -t55 / 0.4e1;
t446 = -t56 / 0.4e1;
t445 = -t57 / 0.4e1;
t441 = mrSges(2,2) * pkin(1);
t440 = mrSges(3,2) * pkin(1);
t53 = pkin(1) * t55;
t54 = pkin(1) * t56;
t439 = pkin(5) * mrSges(2,2);
t52 = t57 * pkin(1);
t438 = -t192 / 0.2e1;
t437 = -t193 / 0.2e1;
t436 = -t194 / 0.2e1;
t435 = -0.3e1 / 0.4e1 * t280;
t434 = m(3) * t280;
t195 = mrSges(2,1) + t265;
t433 = pkin(1) * t195;
t426 = pkin(3) * t280;
t20 = -t53 + t31;
t223 = t251 ^ 2;
t202 = t278 + t280;
t392 = (0.2e1 * t162 + t202) * t226;
t412 = t192 * t81;
t10 = t113 * t192 / (t199 + (t228 * t251 - t377) * pkin(3)) * t392 - (-t20 * t348 + t461 * t31 - (t348 * t454 - t20) * t479 - (-t223 * t451 + t467) * t55) * t55 * t412;
t105 = -g(1) * t183 + g(2) * t186;
t108 = g(1) * t186 + g(2) * t183;
t264 = m(3) * pkin(5);
t327 = -m(3) * qJ(3,3) - mrSges(3,3);
t145 = t264 - t327;
t269 = -pkin(5) / 0.2e1;
t196 = -qJ(3,3) / 0.2e1 + t269;
t262 = Ifges(3,6) / 0.2e1;
t263 = Ifges(3,5) / 0.2e1;
t291 = t80 * t258;
t415 = mrSges(3,2) * t228;
t316 = t181 + t415;
t292 = (mrSges(2,2) + t316) * t476;
t182 = mrSges(3,1) * t228;
t416 = mrSges(3,2) * t227;
t325 = t265 - t416;
t104 = -t182 - t325;
t102 = -mrSges(2,1) + t104;
t149 = t181 + mrSges(2,2);
t128 = t149 + t415;
t266 = m(3) * pkin(1);
t148 = (m(2) * pkin(1)) + mrSges(1,1) + t266;
t305 = t102 * t251 + t128 * t245 - t148;
t339 = t112 * t122 * t473;
t281 = pkin(1) ^ 2;
t374 = t239 * t209;
t290 = (m(2) + m(3)) * t281 + Ifges(2,1) + Ifges(3,2) + Ifges(1,3) - t374 + ((2 * t241 + t480) * pkin(5));
t306 = 0.2e1 * t80;
t323 = pkin(1) * t456 * t228 + (mrSges(2,1) + t325) * t476;
t358 = -0.2e1 * pkin(1) * t149;
t413 = Ifges(3,4) * t227;
t79 = 0.2e1 * t374 + (pkin(2) * t456 + 0.4e1 * t413) * t228 + t434 - 0.2e1 * pkin(2) * t416 + Ifges(2,2) - Ifges(2,1) - t239;
t34 = t79 * t223 + (t306 * t245 + t323) * t251 + (-t245 * t440 - t413) * t448 + t245 * t358 + t233 ^ 2 * m(3) + qJ(3,3) * t455 + t290;
t395 = t112 * t258;
t345 = t55 * t395;
t388 = (-t439 / 0.2e1 + Ifges(2,6) / 0.2e1) * t258;
t397 = pkin(5) * mrSges(2,1) - Ifges(2,5);
t400 = t80 * t223;
t409 = t245 * t79;
t49 = t55 * t433;
t333 = Ifges(2,6) - t439;
t396 = -t195 * pkin(5) + Ifges(2,5);
t58 = (t327 * pkin(2) + t396 - t481) * t245 + t251 * (t142 * t227 + t228 * t460 + t333);
t268 = 0.2e1 * pkin(7);
t169 = cos(t268 + qJ(2,3));
t342 = t73 * t395;
t17 = -t53 - (-t403 / 0.2e1 - t406 / 0.2e1 - t342 / 0.4e1) * t192;
t173 = cos(qJ(2,3) - pkin(7));
t270 = 0.3e1 * qJ(2,3);
t277 = pkin(3) * t278;
t321 = 0.3e1 / 0.4e1 * t280;
t322 = 0.3e1 / 0.4e1 * t278;
t326 = -0.2e1 * t277 - 0.4e1 * t426;
t332 = t112 * t438;
t356 = -0.2e1 * t426;
t357 = -0.2e1 * t432;
t361 = -0.6e1 * t280 - 0.3e1 * t278;
t364 = cos(t210) + t228;
t46 = (t203 ^ 2 + t281) * t55;
t70 = t203 * t395;
t7 = (t169 * t357 + t326 * t172 + t173 * t356 + t251 * t449 + t450) * t332 * t473 + (t295 * t113 * t203 * t474 + ((-t277 * cos(0.3e1 * t211) - t279 * cos(t270)) * t447 - (t373 / 0.2e1 + t370 / 0.2e1) * t70 + (-(-cos(pkin(7) + t270) - t173) * t55 * t321 + (-t31 * t476 + t361 * t447 + t46) * t172) * pkin(3) + ((-t55 * t435 + t46) * t251 - (-cos(t268 + t270) - t169 - 0.2e1 * t251) * t55 * t322 + (-t70 * t163 - 0.2e1 * t17 * t364) * pkin(3) - (t364 * pkin(3) + t251 * t476) * t31) * pkin(2) + (-t31 / 0.2e1 - t17) * (cos(t153) * t278 + cos(t271) * t280 + t202)) * t112) * t192 * t55;
t76 = t104 * t251 + t316 * t245 - t266;
t1 = -t34 * t10 + t58 * t339 - t76 * t7 - 0.4e1 * t345 * t400 - ((pkin(2) * t145 + t397 + t481) * t395 - (t292 + 0.2e1 * t409) * t55) * t251 * t395 - 0.2e1 * (((mrSges(3,2) * t196 + t262) * t395 - mrSges(3,1) * t53) * t228 + ((mrSges(3,1) * t196 + t263) * t395 + mrSges(3,2) * t53) * t227 + t112 * t388 - t49) * t245 * t395 + (-t112 * t291 - t145 * t31) * t454 + (t305 * t105 + t466 * t108) * t252 + t246 * (t105 * t466 - t305 * t108);
t425 = t1 * t192;
t299 = t316 * t251;
t4 = -t76 * t10 - t55 ^ 2 * (-t189 + mrSges(3,3)) + (t105 * t252 - t108 * t246 - t7) * m(3) - 0.2e1 * (-t104 * t245 + t299) * t345;
t424 = t192 * t4;
t106 = -g(1) * t184 + g(2) * t187;
t109 = g(1) * t187 + g(2) * t184;
t21 = -t54 + t32;
t224 = t253 ^ 2;
t411 = t193 * t82;
t11 = t116 * t193 / (t200 + (t228 * t253 - t376) * pkin(3)) * t392 - (-t21 * t347 + t462 * t32 - (t347 * t453 - t21) * t478 - (-t224 * t451 + t467) * t56) * t56 * t411;
t328 = -m(3) * qJ(3,2) - mrSges(3,3);
t146 = t264 - t328;
t197 = -qJ(3,2) / 0.2e1 + t269;
t304 = t102 * t253 + t128 * t247 - t148;
t338 = t115 * t123 * t472;
t394 = t115 * t258;
t344 = t56 * t394;
t35 = t79 * t224 + (t306 * t247 + t323) * t253 + (-t247 * t440 - t413) * t448 + t247 * t358 + t234 ^ 2 * m(3) + qJ(3,2) * t455 + t290;
t399 = t80 * t224;
t408 = t247 * t79;
t50 = t56 * t433;
t59 = (t328 * pkin(2) + t396 - t482) * t247 + t253 * (t143 * t227 + t228 * t459 + t333);
t77 = t104 * t253 + t316 * t247 - t266;
t170 = cos(t268 + qJ(2,2));
t176 = cos(qJ(2,2) - pkin(7));
t341 = t74 * t394;
t18 = -t54 - (-t402 / 0.2e1 - t405 / 0.2e1 - t341 / 0.4e1) * t193;
t272 = 0.3e1 * qJ(2,2);
t331 = t115 * t437;
t363 = cos(t212) + t228;
t47 = (t204 ^ 2 + t281) * t56;
t71 = t204 * t394;
t8 = (t170 * t357 + t326 * t175 + t176 * t356 + t253 * t449 + t450) * t331 * t472 + (t294 * t116 * t204 * t474 + ((-t277 * cos(0.3e1 * t213) - t279 * cos(t272)) * t446 - (t372 / 0.2e1 + t369 / 0.2e1) * t71 + (-(-cos(t272 + pkin(7)) - t176) * t56 * t321 + (-t32 * t476 + t361 * t446 + t47) * t175) * pkin(3) + ((-t56 * t435 + t47) * t253 - (-cos(t272 + t268) - t170 - 0.2e1 * t253) * t56 * t322 + (-t71 * t165 - 0.2e1 * t18 * t363) * pkin(3) - (t363 * pkin(3) + t253 * t476) * t32) * pkin(2) + (-t32 / 0.2e1 - t18) * (cos(t154) * t278 + cos(t273) * t280 + t202)) * t115) * t193 * t56;
t2 = -t35 * t11 + t59 * t338 - t77 * t8 - 0.4e1 * t344 * t399 - ((pkin(2) * t146 + t397 + t482) * t394 - (t292 + 0.2e1 * t408) * t56) * t253 * t394 - 0.2e1 * (((mrSges(3,2) * t197 + t262) * t394 - mrSges(3,1) * t54) * t228 + ((mrSges(3,1) * t197 + t263) * t394 + mrSges(3,2) * t54) * t227 + t115 * t388 - t50) * t247 * t394 + (-t115 * t291 - t146 * t32) * t453 + (t304 * t106 + t465 * t109) * t254 + t248 * (t106 * t465 - t304 * t109);
t423 = t193 * t2;
t298 = t316 * t253;
t5 = -t77 * t11 - t56 ^ 2 * (-t190 + mrSges(3,3)) + (t106 * t254 - t109 * t248 - t8) * m(3) - 0.2e1 * (-t104 * t247 + t298) * t344;
t422 = t193 * t5;
t107 = -g(1) * t185 + g(2) * t188;
t110 = g(1) * t188 + g(2) * t185;
t19 = -t52 + t33;
t225 = t255 ^ 2;
t410 = t194 * t83;
t12 = t119 * t194 / (t201 + (t228 * t255 - t375) * pkin(3)) * t392 - (-t19 * t346 + t463 * t33 - (t346 * t452 - t19) * t477 - (-t225 * t451 + t467) * t57) * t57 * t410;
t329 = -m(3) * qJ(3,1) - mrSges(3,3);
t147 = t264 - t329;
t198 = -qJ(3,1) / 0.2e1 + t269;
t303 = t102 * t255 + t128 * t249 - t148;
t337 = t118 * t124 * t471;
t393 = t118 * t258;
t343 = t57 * t393;
t36 = t79 * t225 + (t306 * t249 + t323) * t255 + (-t249 * t440 - t413) * t448 + t249 * t358 + t235 ^ 2 * m(3) + qJ(3,1) * t455 + t290;
t398 = t80 * t225;
t407 = t249 * t79;
t51 = t57 * t433;
t60 = (t329 * pkin(2) + t396 - t483) * t249 + t255 * (t144 * t227 + t228 * t458 + t333);
t78 = t104 * t255 + t316 * t249 - t266;
t340 = t75 * t393;
t16 = -t52 - (-t401 / 0.2e1 - t404 / 0.2e1 - t340 / 0.4e1) * t194;
t178 = cos(qJ(2,1) + t268);
t180 = cos(qJ(2,1) - pkin(7));
t274 = 0.3e1 * qJ(2,1);
t330 = t118 * t436;
t362 = cos(t214) + t228;
t48 = (t205 ^ 2 + t281) * t57;
t72 = t205 * t393;
t9 = (t178 * t357 + t326 * t179 + t180 * t356 + t255 * t449 + t450) * t330 * t471 + (t293 * t119 * t205 * t474 + ((-t277 * cos(0.3e1 * t215) - t279 * cos(t274)) * t445 - (t371 / 0.2e1 + t368 / 0.2e1) * t72 + (-(-cos(t274 + pkin(7)) - t180) * t57 * t321 + (-t33 * t476 + t361 * t445 + t48) * t179) * pkin(3) + ((-t57 * t435 + t48) * t255 - (-cos(t274 + t268) - t178 - 0.2e1 * t255) * t57 * t322 + (-0.2e1 * t16 * t362 - t72 * t167) * pkin(3) - (t362 * pkin(3) + t255 * t476) * t33) * pkin(2) + (-t33 / 0.2e1 - t16) * (cos(t155) * t278 + cos(t275) * t280 + t202)) * t118) * t194 * t57;
t3 = -t36 * t12 + t60 * t337 - t78 * t9 - 0.4e1 * t343 * t398 - ((pkin(2) * t147 + t397 + t483) * t393 - (t292 + 0.2e1 * t407) * t57) * t255 * t393 - 0.2e1 * (((mrSges(3,2) * t198 + t262) * t393 - mrSges(3,1) * t52) * t228 + ((mrSges(3,1) * t198 + t263) * t393 + mrSges(3,2) * t52) * t227 + t118 * t388 - t51) * t249 * t393 + (-t118 * t291 - t147 * t33) * t452 + (t303 * t107 + t464 * t110) * t256 + t250 * (t107 * t464 - t303 * t110);
t421 = t194 * t3;
t297 = t316 * t255;
t6 = -t78 * t12 - t57 ^ 2 * (-t191 + mrSges(3,3)) + (t107 * t256 - t110 * t250 - t9) * m(3) - 0.2e1 * (-t104 * t249 + t297) * t343;
t420 = t194 * t6;
t242 = xDDP(3);
t387 = t192 * t242;
t243 = xDDP(2);
t386 = t192 * t243;
t244 = xDDP(1);
t385 = t192 * t244;
t384 = t193 * t242;
t383 = t193 * t243;
t382 = t193 * t244;
t381 = t194 * t242;
t380 = t194 * t243;
t379 = t194 * t244;
t367 = -0.2e1 * t265;
t353 = t84 * t412;
t351 = t85 * t411;
t349 = t86 * t410;
t317 = -t416 + t182;
t309 = -t105 * t246 - t108 * t252;
t308 = -t106 * t248 - t109 * t254;
t307 = -t107 * t250 - t110 * t256;
t103 = g(3) * t128;
t101 = g(3) * t102;
t93 = t317 * t475 + Ifges(2,3) + Ifges(3,3) + t434;
t45 = (m(3) * t334 + t78 * t417) * t194;
t44 = (m(3) * t335 + t77 * t418) * t193;
t43 = (m(3) * t336 + t76 * t419) * t192;
t42 = (m(3) * t66 + t78 * t99) * t194;
t41 = (m(3) * t65 + t77 * t98) * t193;
t40 = (m(3) * t64 + t76 * t97) * t192;
t39 = (m(3) * t63 + t78 * t96) * t194;
t38 = (m(3) * t62 + t77 * t95) * t193;
t37 = (m(3) * t61 + t76 * t94) * t192;
t30 = (t36 * t99 + t66 * t78) * t194;
t29 = (t35 * t98 + t65 * t77) * t193;
t28 = (t34 * t97 + t64 * t76) * t192;
t27 = (t36 * t96 + t63 * t78) * t194;
t26 = (t35 * t95 + t62 * t77) * t193;
t25 = (t34 * t94 + t61 * t76) * t192;
t15 = -t36 * t349 + (t75 * t78 * t436 + t60) * t118;
t14 = -t35 * t351 + (t74 * t77 * t437 + t59) * t115;
t13 = -t34 * t353 + (t73 * t76 * t438 + t58) * t112;
t22 = [-(-t27 * t96 - t39 * t63) * t379 - (-t27 * t99 - t39 * t66) * t380 - (-t27 * t417 + (-t39 * t442 + t96 * t60) * t118) * t381 - t96 * t421 - t63 * t420 - (-t26 * t95 - t38 * t62) * t382 - (-t26 * t98 - t38 * t65) * t383 - (-t26 * t418 + (-t38 * t443 + t95 * t59) * t115) * t384 - t95 * t423 - t62 * t422 - (-t25 * t94 - t37 * t61) * t385 - (-t25 * t97 - t37 * t64) * t386 - (-t25 * t419 + (-t37 * t444 + t94 * t58) * t112) * t387 - t94 * t425 - t61 * t424 + (-g(1) + t244) * m(4); -(-t30 * t96 - t42 * t63) * t379 - (-t30 * t99 - t42 * t66) * t380 - (-t30 * t417 + (-t42 * t442 + t99 * t60) * t118) * t381 - t99 * t421 - t66 * t420 - (-t29 * t95 - t41 * t62) * t382 - (-t29 * t98 - t41 * t65) * t383 - (-t29 * t418 + (-t41 * t443 + t98 * t59) * t115) * t384 - t98 * t423 - t65 * t422 - (-t28 * t94 - t40 * t61) * t385 - (-t28 * t97 - t40 * t64) * t386 - (-t28 * t419 + (-t40 * t444 + t97 * t58) * t112) * t387 - t97 * t425 - t64 * t424 + (-g(2) + t243) * m(4); -(t15 * t96 - t45 * t63) * t379 - (t15 * t99 - t45 * t66) * t380 - t3 * t349 + t118 * (-t60 * t12 + t93 * t337 - ((-t33 * t367 - t51) * t249 + (t317 * t249 + t297) * (-t52 - (-t340 - 0.2e1 * t401 - 0.2e1 * t404) * t194) - (-0.2e1 * t398 + (t407 + t441) * t255 + t80) * t57) * t57 + (t307 * t102 + t103) * t249 - (t307 * t128 - t101) * t255) + t75 * t6 * t330 - (t14 * t95 - t44 * t62) * t382 - (t14 * t98 - t44 * t65) * t383 - t2 * t351 + t115 * (-t59 * t11 + t93 * t338 - ((-t32 * t367 - t50) * t247 + (t317 * t247 + t298) * (-t54 - (-t341 - 0.2e1 * t402 - 0.2e1 * t405) * t193) - (-0.2e1 * t399 + (t408 + t441) * t253 + t80) * t56) * t56 + (t308 * t102 + t103) * t247 - (t308 * t128 - t101) * t253) + t74 * t5 * t331 - (t13 * t94 - t43 * t61) * t385 - (t13 * t97 - t43 * t64) * t386 - t1 * t353 + t112 * (-t58 * t10 + t93 * t339 - ((-t31 * t367 - t49) * t245 + (t317 * t245 + t299) * (-t53 - (-t342 - 0.2e1 * t403 - 0.2e1 * t406) * t192) - (-0.2e1 * t400 + (t409 + t441) * t251 + t80) * t55) * t55 + (t309 * t102 + t103) * t245 - (t309 * t128 - t101) * t251) + t73 * t4 * t332 - g(3) * m(4) + (-t15 * t349 + (t118 * t93 - (t60 * t417 - t45 * t442) * t194) * t118 - t14 * t351 + (t115 * t93 - (t59 * t418 - t44 * t443) * t193) * t115 - t13 * t353 + (t112 * t93 - (t58 * t419 - t43 * t444) * t192) * t112 + m(4)) * t242;];
tauX  = t22;
