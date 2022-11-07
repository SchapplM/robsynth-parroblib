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
% Datum: 2022-11-07 13:12
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
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
% StartTime: 2022-11-07 13:11:07
% EndTime: 2022-11-07 13:11:15
% DurationCPUTime: 8.25s
% Computational Cost: add. (29634->565), mult. (33291->883), div. (5214->9), fcn. (23079->62), ass. (0->362)
t219 = qJ(3,1) + pkin(5);
t128 = mrSges(3,2) * t219 - Ifges(3,6);
t212 = sin(pkin(7));
t213 = cos(pkin(7));
t418 = mrSges(3,1) * t219 - Ifges(3,5);
t431 = -t128 * t212 + t418 * t213;
t218 = qJ(3,2) + pkin(5);
t127 = mrSges(3,2) * t218 - Ifges(3,6);
t419 = mrSges(3,1) * t218 - Ifges(3,5);
t430 = -t127 * t212 + t419 * t213;
t217 = qJ(3,3) + pkin(5);
t126 = mrSges(3,2) * t217 - Ifges(3,6);
t420 = mrSges(3,1) * t217 - Ifges(3,5);
t429 = -t126 * t212 + t420 * t213;
t428 = (m(2) * pkin(5));
t427 = 2 * pkin(1);
t190 = -pkin(6) - t219;
t179 = 0.1e1 / t190;
t236 = cos(qJ(2,1));
t186 = t236 * pkin(2);
t200 = qJ(2,1) + pkin(7);
t164 = cos(t200);
t391 = pkin(3) * t164;
t111 = t186 + t391;
t106 = 0.1e1 / t111;
t239 = xDP(3);
t357 = t106 * t239;
t230 = sin(qJ(2,1));
t114 = pkin(2) * t230 + pkin(3) * sin(t200);
t256 = 0.2e1 * qJ(2,1);
t199 = t256 + pkin(7);
t152 = sin(t199);
t408 = pkin(2) * pkin(3);
t336 = 0.2e1 * t408;
t261 = pkin(2) ^ 2;
t339 = t261 * sin(t256);
t141 = 0.2e1 * t200;
t259 = pkin(3) ^ 2;
t342 = t259 * sin(t141);
t277 = t152 * t336 + t339 + t342;
t72 = t114 * t427 + t277;
t319 = t72 * t357;
t241 = xDP(1);
t216 = legFrame(1,3);
t170 = sin(t216);
t173 = cos(t216);
t144 = t186 + pkin(1);
t231 = sin(qJ(1,1));
t237 = cos(qJ(1,1));
t80 = t144 * t231 + t190 * t237;
t83 = t144 * t237 - t190 * t231;
t90 = -t170 * t231 + t173 * t237;
t66 = -t170 * t80 + t173 * t83 + t90 * t391;
t371 = t66 * t241;
t240 = xDP(2);
t89 = t170 * t237 + t173 * t231;
t65 = t170 * t83 + t173 * t80 + t89 * t391;
t372 = t65 * t240;
t33 = (t371 + t372 + t319 / 0.2e1) * t179;
t189 = -pkin(6) - t218;
t178 = 0.1e1 / t189;
t234 = cos(qJ(2,2));
t185 = t234 * pkin(2);
t198 = qJ(2,2) + pkin(7);
t160 = cos(t198);
t392 = pkin(3) * t160;
t110 = t185 + t392;
t104 = 0.1e1 / t110;
t360 = t104 * t239;
t228 = sin(qJ(2,2));
t113 = pkin(2) * t228 + pkin(3) * sin(t198);
t254 = 0.2e1 * qJ(2,2);
t197 = t254 + pkin(7);
t150 = sin(t197);
t340 = t261 * sin(t254);
t140 = 0.2e1 * t198;
t343 = t259 * sin(t140);
t278 = t150 * t336 + t340 + t343;
t71 = t113 * t427 + t278;
t320 = t71 * t360;
t215 = legFrame(2,3);
t169 = sin(t215);
t172 = cos(t215);
t143 = t185 + pkin(1);
t229 = sin(qJ(1,2));
t235 = cos(qJ(1,2));
t79 = t143 * t229 + t189 * t235;
t82 = t143 * t235 - t189 * t229;
t88 = -t169 * t229 + t172 * t235;
t64 = -t169 * t79 + t172 * t82 + t88 * t392;
t373 = t64 * t241;
t87 = t169 * t235 + t172 * t229;
t63 = t169 * t82 + t172 * t79 + t87 * t392;
t374 = t63 * t240;
t32 = (t373 + t374 + t320 / 0.2e1) * t178;
t188 = -pkin(6) - t217;
t177 = 0.1e1 / t188;
t232 = cos(qJ(2,3));
t184 = t232 * pkin(2);
t196 = qJ(2,3) + pkin(7);
t156 = cos(t196);
t393 = pkin(3) * t156;
t109 = t184 + t393;
t102 = 0.1e1 / t109;
t363 = t102 * t239;
t226 = sin(qJ(2,3));
t112 = pkin(2) * t226 + pkin(3) * sin(t196);
t252 = 0.2e1 * qJ(2,3);
t195 = t252 + pkin(7);
t148 = sin(t195);
t341 = t261 * sin(t252);
t139 = 0.2e1 * t196;
t344 = t259 * sin(t139);
t279 = t148 * t336 + t341 + t344;
t70 = t112 * t427 + t279;
t321 = t70 * t363;
t214 = legFrame(3,3);
t168 = sin(t214);
t171 = cos(t214);
t142 = t184 + pkin(1);
t227 = sin(qJ(1,3));
t233 = cos(qJ(1,3));
t78 = t142 * t227 + t188 * t233;
t81 = t142 * t233 - t188 * t227;
t86 = -t168 * t227 + t171 * t233;
t62 = -t168 * t78 + t171 * t81 + t86 * t393;
t375 = t62 * t241;
t85 = t168 * t233 + t171 * t227;
t61 = t168 * t81 + t171 * t78 + t85 * t393;
t376 = t61 * t240;
t31 = (t375 + t376 + t321 / 0.2e1) * t177;
t166 = t212 * mrSges(3,1);
t194 = t213 ^ 2;
t220 = Ifges(3,2) - Ifges(3,1);
t77 = (-pkin(2) * mrSges(3,2) - t212 * t220) * t213 - pkin(2) * t166 + 0.2e1 * Ifges(3,4) * t194 + Ifges(2,4) - Ifges(3,4);
t174 = (t217 * m(3));
t222 = (mrSges(2,3) + mrSges(3,3));
t302 = -mrSges(1,2) + t222 + t428;
t423 = -t302 - t174;
t175 = (t218 * m(3));
t422 = -t302 - t175;
t176 = (t219 * m(3));
t421 = -t302 - t176;
t416 = 0.2e1 * mrSges(3,1);
t415 = 2 * mrSges(3,3);
t211 = t239 ^ 2;
t414 = 0.1e1 / t109 ^ 2;
t413 = 0.1e1 / t110 ^ 2;
t412 = 0.1e1 / t111 ^ 2;
t326 = t213 * t408;
t411 = -0.4e1 * pkin(1) * (t326 + t261 / 0.2e1 + t259 / 0.2e1);
t260 = pkin(2) * t261;
t394 = pkin(2) * t259;
t410 = -0.2e1 * t260 - 0.4e1 * t394;
t409 = 0.2e1 * t213;
t246 = m(3) * pkin(2);
t356 = t112 * t102;
t55 = (t239 * t356 + t240 * t85 + t241 * t86) * t177;
t407 = -t55 / 0.4e1;
t355 = t113 * t104;
t56 = (t239 * t355 + t240 * t87 + t241 * t88) * t178;
t406 = -t56 / 0.4e1;
t354 = t114 * t106;
t57 = (t239 * t354 + t240 * t89 + t241 * t90) * t179;
t405 = -t57 / 0.4e1;
t404 = t70 / 0.2e1;
t403 = t71 / 0.2e1;
t402 = t72 / 0.2e1;
t401 = mrSges(2,2) * pkin(1);
t400 = mrSges(3,2) * pkin(1);
t53 = pkin(1) * t55;
t54 = pkin(1) * t56;
t399 = pkin(5) * mrSges(2,2);
t52 = t57 * pkin(1);
t398 = -t239 / 0.2e1;
t397 = -0.3e1 / 0.4e1 * t261;
t396 = m(3) * t261;
t180 = mrSges(2,1) + t246;
t395 = pkin(1) * t180;
t390 = pkin(3) * t261;
t103 = t102 * t414;
t155 = cos(t195);
t187 = t259 + t261;
t286 = cos(t139) * t259 + cos(t252) * t261 + t187;
t312 = -0.2e1 * t326;
t337 = -0.2e1 * t408;
t353 = (0.2e1 * t326 + t187) * t211;
t365 = t102 * t177;
t10 = t103 * t177 * t353 + (-(t155 * t337 - t286 + t312) * t55 / 0.2e1 + (-0.2e1 * t31 + t53) * t109) * t55 * t365;
t245 = m(3) * pkin(5);
t305 = -m(3) * qJ(3,3) - mrSges(3,3);
t132 = t245 - t305;
t250 = -pkin(5) / 0.2e1;
t181 = -qJ(3,3) / 0.2e1 + t250;
t243 = Ifges(3,6) / 0.2e1;
t244 = Ifges(3,5) / 0.2e1;
t272 = t77 * t239;
t382 = mrSges(3,2) * t213;
t294 = t166 + t382;
t276 = (mrSges(2,2) + t294) * t427;
t136 = t166 + mrSges(2,2);
t115 = t136 + t382;
t247 = m(3) * pkin(1);
t135 = (m(2) * pkin(1)) + mrSges(1,1) + t247;
t167 = mrSges(3,1) * t213;
t383 = mrSges(3,2) * t212;
t303 = t246 - t383;
t94 = -t167 - t303;
t92 = -mrSges(2,1) + t94;
t289 = t115 * t226 + t232 * t92 - t135;
t317 = t103 * t211 * t112;
t324 = t55 * t363;
t208 = t232 ^ 2;
t262 = pkin(1) ^ 2;
t345 = t220 * t194;
t271 = (m(2) + m(3)) * t262 + Ifges(2,1) + Ifges(3,2) + Ifges(1,3) - t345 + ((2 * t222 + t428) * pkin(5));
t290 = 0.2e1 * t77;
t301 = pkin(1) * t416 * t213 + (mrSges(2,1) + t303) * t427;
t329 = -0.2e1 * pkin(1) * t136;
t380 = Ifges(3,4) * t212;
t76 = 0.2e1 * t345 + (pkin(2) * t416 + 0.4e1 * t380) * t213 + t396 - 0.2e1 * pkin(2) * t383 + Ifges(2,2) - Ifges(2,1) - t220;
t34 = t76 * t208 + (t226 * t290 + t301) * t232 + (-t226 * t400 - t380) * t409 + t226 * t329 + t217 ^ 2 * m(3) + qJ(3,3) * t415 + t271;
t352 = (-t399 / 0.2e1 + Ifges(2,6) / 0.2e1) * t239;
t367 = pkin(5) * mrSges(2,1) - Ifges(2,5);
t370 = t77 * t208;
t379 = t226 * t76;
t49 = t55 * t395;
t311 = Ifges(2,6) - t399;
t366 = -t180 * pkin(5) + Ifges(2,5);
t58 = (t305 * pkin(2) + t366 - t429) * t226 + t232 * (-t126 * t213 - t212 * t420 + t311);
t249 = 0.2e1 * pkin(7);
t154 = cos(t249 + qJ(2,3));
t157 = cos(qJ(2,3) - pkin(7));
t17 = -t53 - (-t375 / 0.2e1 - t376 / 0.2e1 - t321 / 0.4e1) * t177;
t251 = 0.3e1 * qJ(2,3);
t258 = pkin(3) * t259;
t299 = 0.3e1 / 0.4e1 * t261;
t300 = 0.3e1 / 0.4e1 * t259;
t304 = -0.2e1 * t258 - 0.4e1 * t390;
t310 = -t365 / 0.2e1;
t327 = -0.2e1 * t390;
t328 = -0.2e1 * t394;
t332 = -0.6e1 * t261 - 0.3e1 * t259;
t335 = t155 + t213;
t46 = (t188 ^ 2 + t262) * t55;
t67 = t188 * t363;
t7 = (t154 * t328 + t304 * t156 + t157 * t327 + t232 * t410 + t411) * t211 * t414 * t310 + (t279 * t414 * t188 * t398 + ((-t258 * cos(0.3e1 * t196) - t260 * cos(t251)) * t407 - (t344 / 0.2e1 + t341 / 0.2e1) * t67 + (-(-cos(pkin(7) + t251) - t157) * t55 * t299 + (-t31 * t427 + t332 * t407 + t46) * t156) * pkin(3) + ((-t55 * t397 + t46) * t232 - (-cos(t249 + t251) - t154 - 0.2e1 * t232) * t55 * t300 + (-t67 * t148 - 0.2e1 * t17 * t335) * pkin(3) - (t335 * pkin(3) + t232 * t427) * t31) * pkin(2) + (-t31 / 0.2e1 - t17) * t286) * t102) * t177 * t55;
t73 = t226 * t294 + t232 * t94 - t247;
t95 = -g(1) * t168 + g(2) * t171;
t98 = g(1) * t171 + g(2) * t168;
t389 = (-t34 * t10 + t58 * t317 - t73 * t7 - 0.4e1 * t324 * t370 - ((pkin(2) * t132 + t367 + t429) * t363 - (t276 + 0.2e1 * t379) * t55) * t232 * t363 - 0.2e1 * (((mrSges(3,2) * t181 + t243) * t363 - mrSges(3,1) * t53) * t213 + ((mrSges(3,1) * t181 + t244) * t363 + mrSges(3,2) * t53) * t212 + t102 * t352 - t49) * t226 * t363 - 0.2e1 * t55 * (-t102 * t272 - t132 * t31) + (t289 * t95 + t423 * t98) * t233 + (-t289 * t98 + t95 * t423) * t227) * t177;
t283 = t294 * t232;
t4 = -t73 * t10 - (t174 + mrSges(3,3)) * t55 ^ 2 + (-t227 * t98 + t233 * t95 - t7) * m(3) - 0.2e1 * (-t226 * t94 + t283) * t324;
t388 = t177 * t4;
t105 = t104 * t413;
t158 = cos(t197);
t285 = cos(t140) * t259 + cos(t254) * t261 + t187;
t362 = t104 * t178;
t11 = t105 * t178 * t353 + (-(t158 * t337 - t285 + t312) * t56 / 0.2e1 + (-0.2e1 * t32 + t54) * t110) * t56 * t362;
t306 = -m(3) * qJ(3,2) - mrSges(3,3);
t133 = t245 - t306;
t182 = -qJ(3,2) / 0.2e1 + t250;
t288 = t115 * t228 + t234 * t92 - t135;
t315 = t105 * t211 * t113;
t323 = t56 * t360;
t209 = t234 ^ 2;
t35 = t76 * t209 + (t228 * t290 + t301) * t234 + (-t228 * t400 - t380) * t409 + t228 * t329 + t218 ^ 2 * m(3) + qJ(3,2) * t415 + t271;
t369 = t77 * t209;
t378 = t228 * t76;
t50 = t56 * t395;
t59 = (t306 * pkin(2) + t366 - t430) * t228 + t234 * (-t127 * t213 - t212 * t419 + t311);
t74 = t228 * t294 + t234 * t94 - t247;
t159 = cos(qJ(2,2) + t249);
t161 = cos(qJ(2,2) - pkin(7));
t18 = -t54 - (-t373 / 0.2e1 - t374 / 0.2e1 - t320 / 0.4e1) * t178;
t253 = 0.3e1 * qJ(2,2);
t309 = -t362 / 0.2e1;
t334 = t158 + t213;
t47 = (t189 ^ 2 + t262) * t56;
t68 = t189 * t360;
t8 = (t159 * t328 + t304 * t160 + t161 * t327 + t234 * t410 + t411) * t211 * t413 * t309 + (t278 * t413 * t189 * t398 + ((-t258 * cos(0.3e1 * t198) - t260 * cos(t253)) * t406 - (t343 / 0.2e1 + t340 / 0.2e1) * t68 + (-(-cos(t253 + pkin(7)) - t161) * t56 * t299 + (-t32 * t427 + t332 * t406 + t47) * t160) * pkin(3) + ((-t56 * t397 + t47) * t234 - (-cos(t253 + t249) - t159 - 0.2e1 * t234) * t56 * t300 + (-t68 * t150 - 0.2e1 * t18 * t334) * pkin(3) - (t334 * pkin(3) + t234 * t427) * t32) * pkin(2) + (-t32 / 0.2e1 - t18) * t285) * t104) * t178 * t56;
t96 = -g(1) * t169 + g(2) * t172;
t99 = g(1) * t172 + g(2) * t169;
t387 = t178 * (-t35 * t11 + t59 * t315 - t74 * t8 - 0.4e1 * t323 * t369 - ((pkin(2) * t133 + t367 + t430) * t360 - (t276 + 0.2e1 * t378) * t56) * t234 * t360 - 0.2e1 * (((mrSges(3,2) * t182 + t243) * t360 - mrSges(3,1) * t54) * t213 + ((mrSges(3,1) * t182 + t244) * t360 + mrSges(3,2) * t54) * t212 + t104 * t352 - t50) * t228 * t360 - 0.2e1 * t56 * (-t104 * t272 - t133 * t32) + (t288 * t96 + t422 * t99) * t235 + (-t288 * t99 + t96 * t422) * t229);
t282 = t294 * t234;
t5 = -t74 * t11 - (t175 + mrSges(3,3)) * t56 ^ 2 + (-t229 * t99 + t235 * t96 - t8) * m(3) - 0.2e1 * (-t228 * t94 + t282) * t323;
t386 = t178 * t5;
t100 = g(1) * t173 + g(2) * t170;
t107 = t106 * t412;
t162 = cos(t199);
t284 = cos(t141) * t259 + cos(t256) * t261 + t187;
t359 = t106 * t179;
t12 = t107 * t179 * t353 + (-(t162 * t337 - t284 + t312) * t57 / 0.2e1 + (-0.2e1 * t33 + t52) * t111) * t57 * t359;
t307 = -m(3) * qJ(3,1) - mrSges(3,3);
t134 = t245 - t307;
t183 = -qJ(3,1) / 0.2e1 + t250;
t287 = t115 * t230 + t236 * t92 - t135;
t313 = t107 * t211 * t114;
t322 = t57 * t357;
t210 = t236 ^ 2;
t36 = t76 * t210 + (t230 * t290 + t301) * t236 + (-t230 * t400 - t380) * t409 + t230 * t329 + t219 ^ 2 * m(3) + qJ(3,1) * t415 + t271;
t368 = t77 * t210;
t377 = t230 * t76;
t51 = t57 * t395;
t60 = (t307 * pkin(2) + t366 - t431) * t230 + t236 * (-t128 * t213 - t212 * t418 + t311);
t75 = t230 * t294 + t236 * t94 - t247;
t16 = -t52 - (-t371 / 0.2e1 - t372 / 0.2e1 - t319 / 0.4e1) * t179;
t163 = cos(qJ(2,1) + t249);
t165 = cos(qJ(2,1) - pkin(7));
t255 = 0.3e1 * qJ(2,1);
t308 = -t359 / 0.2e1;
t333 = t162 + t213;
t48 = (t190 ^ 2 + t262) * t57;
t69 = t190 * t357;
t9 = (t163 * t328 + t304 * t164 + t165 * t327 + t236 * t410 + t411) * t211 * t412 * t308 + (t277 * t412 * t190 * t398 + ((-t258 * cos(0.3e1 * t200) - t260 * cos(t255)) * t405 - (t342 / 0.2e1 + t339 / 0.2e1) * t69 + (-(-cos(pkin(7) + t255) - t165) * t57 * t299 + (-t33 * t427 + t332 * t405 + t48) * t164) * pkin(3) + ((-t57 * t397 + t48) * t236 - (-cos(t255 + t249) - t163 - 0.2e1 * t236) * t57 * t300 + (-t69 * t152 - 0.2e1 * t16 * t333) * pkin(3) - (t333 * pkin(3) + t236 * t427) * t33) * pkin(2) + (-t33 / 0.2e1 - t16) * t284) * t106) * t179 * t57;
t97 = -g(1) * t170 + g(2) * t173;
t385 = t179 * (-t36 * t12 + t60 * t313 - t75 * t9 - 0.4e1 * t322 * t368 - ((pkin(2) * t134 + t367 + t431) * t357 - (t276 + 0.2e1 * t377) * t57) * t236 * t357 - 0.2e1 * (((mrSges(3,2) * t183 + t243) * t357 - mrSges(3,1) * t52) * t213 + ((mrSges(3,1) * t183 + t244) * t357 + mrSges(3,2) * t52) * t212 + t106 * t352 - t51) * t230 * t357 - 0.2e1 * t57 * (-t106 * t272 - t134 * t33) + (t421 * t100 + t287 * t97) * t237 + (-t287 * t100 + t97 * t421) * t231);
t281 = t294 * t236;
t6 = -t75 * t12 - (t176 + mrSges(3,3)) * t57 ^ 2 + (-t100 * t231 + t237 * t97 - t9) * m(3) - 0.2e1 * (-t230 * t94 + t281) * t322;
t384 = t179 * t6;
t223 = xDDP(3);
t364 = t102 * t223;
t361 = t104 * t223;
t358 = t106 * t223;
t224 = xDDP(2);
t351 = t177 * t224;
t225 = xDDP(1);
t350 = t177 * t225;
t349 = t178 * t224;
t348 = t178 * t225;
t347 = t179 * t224;
t346 = t179 * t225;
t338 = -0.2e1 * t246;
t318 = t177 * t364;
t316 = t178 * t361;
t314 = t179 * t358;
t295 = -t383 + t167;
t293 = -t227 * t95 - t233 * t98;
t292 = -t229 * t96 - t235 * t99;
t291 = -t100 * t237 - t231 * t97;
t93 = g(3) * t115;
t91 = g(3) * t92;
t84 = 0.2e1 * pkin(2) * t295 + Ifges(2,3) + Ifges(3,3) + t396;
t45 = (m(3) * t402 + t114 * t75) * t359;
t44 = (m(3) * t403 + t113 * t74) * t362;
t43 = (m(3) * t404 + t112 * t73) * t365;
t42 = (m(3) * t66 + t75 * t90) * t179;
t41 = (m(3) * t65 + t75 * t89) * t179;
t40 = (m(3) * t64 + t74 * t88) * t178;
t39 = (m(3) * t63 + t74 * t87) * t178;
t38 = (m(3) * t62 + t73 * t86) * t177;
t37 = (m(3) * t61 + t73 * t85) * t177;
t30 = (t36 * t90 + t66 * t75) * t179;
t29 = (t36 * t89 + t65 * t75) * t179;
t28 = (t35 * t88 + t64 * t74) * t178;
t27 = (t35 * t87 + t63 * t74) * t178;
t26 = (t34 * t86 + t62 * t73) * t177;
t25 = (t34 * t85 + t61 * t73) * t177;
t15 = (t60 - (t114 * t36 + t75 * t402) * t179) * t106;
t14 = (t59 - (t113 * t35 + t74 * t403) * t178) * t104;
t13 = (t58 - (t112 * t34 + t73 * t404) * t177) * t102;
t1 = [-(-t30 * t90 - t42 * t66) * t346 - (-t30 * t89 - t42 * t65) * t347 - (-t30 * t114 - t42 * t402 + t90 * t60) * t314 - t90 * t385 - t66 * t384 - (-t28 * t88 - t40 * t64) * t348 - (-t28 * t87 - t40 * t63) * t349 - (-t28 * t113 - t40 * t403 + t88 * t59) * t316 - t88 * t387 - t64 * t386 - (-t26 * t86 - t38 * t62) * t350 - (-t26 * t85 - t38 * t61) * t351 - (-t26 * t112 - t38 * t404 + t86 * t58) * t318 - t86 * t389 - t62 * t388 + (-g(1) + t225) * m(4); -(-t29 * t90 - t41 * t66) * t346 - (-t29 * t89 - t41 * t65) * t347 - (-t29 * t114 - t41 * t402 + t89 * t60) * t314 - t89 * t385 - t65 * t384 - (-t27 * t88 - t39 * t64) * t348 - (-t27 * t87 - t39 * t63) * t349 - (-t27 * t113 - t39 * t403 + t87 * t59) * t316 - t87 * t387 - t63 * t386 - (-t25 * t86 - t37 * t62) * t350 - (-t25 * t85 - t37 * t61) * t351 - (-t25 * t112 - t37 * t404 + t85 * t58) * t318 - t85 * t389 - t61 * t388 + (-g(2) + t224) * m(4); -(t15 * t90 - t45 * t66) * t346 - (t15 * t89 - t45 * t65) * t347 + (t106 * t84 - (-t45 * t402 + (t106 * t60 + t15) * t114) * t179) * t358 - t354 * t385 + t106 * (-t60 * t12 + t84 * t313 - ((-t33 * t338 - t51) * t230 + (t230 * t295 + t281) * (-t52 - (-t319 - 0.2e1 * t371 - 0.2e1 * t372) * t179) - (-0.2e1 * t368 + (t377 + t401) * t236 + t77) * t57) * t57 + (t291 * t92 + t93) * t230 - t236 * (t115 * t291 - t91)) + t72 * t6 * t308 - (t14 * t88 - t44 * t64) * t348 - (t14 * t87 - t44 * t63) * t349 + (t104 * t84 - (-t44 * t403 + (t104 * t59 + t14) * t113) * t178) * t361 - t355 * t387 + t104 * (-t59 * t11 + t84 * t315 - ((-t32 * t338 - t50) * t228 + (t228 * t295 + t282) * (-t54 - (-t320 - 0.2e1 * t373 - 0.2e1 * t374) * t178) - (-0.2e1 * t369 + (t378 + t401) * t234 + t77) * t56) * t56 + (t292 * t92 + t93) * t228 - t234 * (t115 * t292 - t91)) + t71 * t5 * t309 - (t13 * t86 - t43 * t62) * t350 - (t13 * t85 - t43 * t61) * t351 + (t102 * t84 - (-t43 * t404 + (t102 * t58 + t13) * t112) * t177) * t364 - t356 * t389 + t102 * (-t58 * t10 + t84 * t317 - ((-t31 * t338 - t49) * t226 + (t226 * t295 + t283) * (-t53 - (-t321 - 0.2e1 * t375 - 0.2e1 * t376) * t177) - (-0.2e1 * t370 + (t379 + t401) * t232 + t77) * t55) * t55 + (t293 * t92 + t93) * t226 - t232 * (t115 * t293 - t91)) + t70 * t4 * t310 + (-g(3) + t223) * m(4);];
tauX  = t1;
