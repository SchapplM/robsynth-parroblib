% Calculate vector of inverse dynamics forces for parallel robot
% P3RPRRR9V1G1A0
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
% Datum: 2020-08-06 18:48
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRRR9V1G1A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G1A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:46:57
% EndTime: 2020-08-06 18:47:06
% DurationCPUTime: 7.70s
% Computational Cost: add. (20409->574), mult. (18900->902), div. (4026->11), fcn. (15345->53), ass. (0->361)
t211 = 2 * pkin(7);
t159 = t211 + qJ(3,1);
t194 = cos(qJ(3,1));
t415 = (cos(t159) + t194) * pkin(2);
t158 = t211 + qJ(3,2);
t192 = cos(qJ(3,2));
t414 = (cos(t158) + t192) * pkin(2);
t157 = t211 + qJ(3,3);
t190 = cos(qJ(3,3));
t413 = (cos(t157) + t190) * pkin(2);
t162 = pkin(7) + qJ(3,1);
t113 = 2 * t162;
t412 = (cos(t113) + 0.1e1) * pkin(3);
t161 = pkin(7) + qJ(3,2);
t112 = 2 * t161;
t411 = (cos(t112) + 0.1e1) * pkin(3);
t160 = pkin(7) + qJ(3,3);
t111 = 2 * t160;
t410 = (cos(t111) + 0.1e1) * pkin(3);
t120 = sin(t160);
t184 = sin(qJ(3,3));
t390 = 2 * pkin(1);
t97 = sin(t111);
t64 = t120 * t390 + pkin(3) * t97 + (sin(t157) + t184) * pkin(2);
t383 = t64 / 0.2e1;
t121 = sin(t161);
t186 = sin(qJ(3,2));
t98 = sin(t112);
t65 = t121 * t390 + pkin(3) * t98 + (sin(t158) + t186) * pkin(2);
t382 = t65 / 0.2e1;
t122 = sin(t162);
t188 = sin(qJ(3,1));
t99 = sin(t113);
t66 = t122 * t390 + pkin(3) * t99 + (sin(t159) + t188) * pkin(2);
t381 = t66 / 0.2e1;
t409 = t415 / 0.2e1 + t412 / 0.2e1;
t408 = t414 / 0.2e1 + t411 / 0.2e1;
t407 = t413 / 0.2e1 + t410 / 0.2e1;
t403 = (m(2) * qJ(2,1));
t402 = (m(2) * qJ(2,2));
t401 = (m(2) * qJ(2,3));
t200 = (m(2) + m(3));
t149 = pkin(1) * t200;
t176 = (pkin(5) + qJ(2,3));
t150 = -pkin(6) - t176;
t400 = t150 * t97;
t177 = (pkin(5) + qJ(2,2));
t151 = -pkin(6) - t177;
t399 = t151 * t98;
t178 = (pkin(5) + qJ(2,1));
t152 = -pkin(6) - t178;
t398 = t152 * t99;
t169 = t194 ^ 2;
t352 = t169 * Ifges(3,4);
t397 = 0.2e1 * t352 - Ifges(3,4);
t168 = t192 ^ 2;
t353 = t168 * Ifges(3,4);
t396 = 0.2e1 * t353 - Ifges(3,4);
t167 = t190 ^ 2;
t354 = t167 * Ifges(3,4);
t395 = 0.2e1 * t354 - Ifges(3,4);
t197 = xDP(3);
t170 = t197 ^ 2;
t218 = 0.1e1 / pkin(3);
t394 = t170 * t218;
t106 = mrSges(3,1) * t176 - Ifges(3,5);
t107 = mrSges(3,1) * t177 - Ifges(3,5);
t108 = mrSges(3,1) * t178 - Ifges(3,5);
t103 = mrSges(3,2) * t176 - Ifges(3,6);
t104 = mrSges(3,2) * t177 - Ifges(3,6);
t105 = mrSges(3,2) * t178 - Ifges(3,6);
t126 = cos(t160);
t114 = 0.1e1 / t126;
t127 = cos(t161);
t116 = 0.1e1 / t127;
t128 = cos(t162);
t118 = 0.1e1 / t128;
t389 = 0.1e1 / t126 ^ 2;
t388 = 0.1e1 / t127 ^ 2;
t387 = 0.1e1 / t128 ^ 2;
t172 = cos(pkin(7));
t156 = t172 ^ 2;
t386 = -0.4e1 * t156;
t171 = sin(pkin(7));
t385 = -0.2e1 * t171;
t180 = (mrSges(2,3) + mrSges(3,3));
t384 = 2 * t180;
t380 = mrSges(3,1) * g(3);
t379 = pkin(1) * mrSges(3,2);
t142 = 0.1e1 / t150;
t198 = xDP(2);
t199 = xDP(1);
t334 = t114 * t197;
t173 = legFrame(3,3);
t130 = sin(t173);
t133 = cos(t173);
t185 = sin(qJ(1,3));
t191 = cos(qJ(1,3));
t73 = -t130 * t185 + t133 * t191;
t76 = t130 * t191 + t133 * t185;
t55 = (t120 * t334 + t198 * t76 + t199 * t73) * t142;
t53 = pkin(1) * t55;
t143 = 0.1e1 / t151;
t331 = t116 * t197;
t174 = legFrame(2,3);
t131 = sin(t174);
t134 = cos(t174);
t187 = sin(qJ(1,2));
t193 = cos(qJ(1,2));
t74 = -t131 * t187 + t134 * t193;
t77 = t131 * t193 + t134 * t187;
t56 = (t121 * t331 + t198 * t77 + t199 * t74) * t143;
t54 = pkin(1) * t56;
t378 = pkin(2) * mrSges(3,1);
t377 = pkin(2) * mrSges(3,2);
t144 = 0.1e1 / t152;
t328 = t118 * t197;
t175 = legFrame(1,3);
t132 = sin(t175);
t135 = cos(t175);
t189 = sin(qJ(1,1));
t195 = cos(qJ(1,1));
t75 = -t132 * t189 + t135 * t195;
t78 = t132 * t195 + t135 * t189;
t57 = (t122 * t328 + t198 * t78 + t199 * t75) * t144;
t52 = t57 * pkin(1);
t376 = -t197 / 0.2e1;
t375 = t200 / 0.2e1;
t275 = t66 * t328;
t371 = pkin(3) * t128;
t129 = t172 * pkin(2);
t110 = t129 + pkin(1);
t69 = t110 * t189 + t152 * t195;
t72 = t110 * t195 - t152 * t189;
t51 = t132 * t72 + t135 * t69 + t78 * t371;
t338 = t51 * t198;
t48 = -t132 * t69 + t135 * t72 + t75 * t371;
t341 = t48 * t199;
t228 = t275 + 0.2e1 * t338 + 0.2e1 * t341;
t25 = t228 * t144 - t52;
t374 = mrSges(3,2) * t25;
t373 = pkin(3) * t126;
t372 = pkin(3) * t127;
t115 = t114 * t389;
t277 = t64 * t334;
t67 = t110 * t185 + t150 * t191;
t70 = t110 * t191 - t150 * t185;
t49 = t130 * t70 + t133 * t67 + t76 * t373;
t340 = t49 * t198;
t46 = -t130 * t67 + t133 * t70 + t73 * t373;
t343 = t46 * t199;
t230 = t277 + 0.2e1 * t340 + 0.2e1 * t343;
t336 = t114 * t142;
t265 = -t336 / 0.2e1;
t28 = (t343 + t340 + t277 / 0.2e1) * t142;
t10 = (-t115 * t394 - (-t28 + t114 * t55 * t407 + (t53 * t114 + t230 * t265) * t126) * t55) * t142;
t147 = -t377 / 0.2e1;
t148 = t378 / 0.4e1;
t201 = Ifges(3,6) / 0.2e1;
t202 = Ifges(3,5) / 0.4e1;
t212 = -pkin(5) / 0.4e1;
t213 = -pkin(5) / 0.2e1;
t285 = -mrSges(1,2) + t180;
t235 = m(3) * t176 + t285 + t401;
t109 = mrSges(1,1) + t149;
t361 = mrSges(3,2) * t184;
t249 = mrSges(3,1) * t190 - t361;
t136 = t184 * mrSges(3,1);
t303 = mrSges(3,2) * t190 + t136;
t337 = m(3) * pkin(2) + mrSges(2,1);
t246 = t171 * (mrSges(2,2) + t303) + t172 * (-t249 - t337);
t240 = -t109 + t246;
t317 = t170 / pkin(3) ^ 2;
t253 = t115 * t120 * t317;
t179 = Ifges(3,1) - Ifges(3,2);
t315 = t179 * t167;
t351 = t171 * t55;
t258 = t315 * t351;
t306 = t218 * t197;
t261 = 0.4e1 * t172 * t306;
t269 = t171 * t306;
t273 = t114 * t306;
t280 = -t378 / 0.2e1;
t281 = -t379 / 0.2e1;
t288 = pkin(2) * t361;
t284 = (t288 - t179) * t351;
t298 = mrSges(3,1) * t53;
t300 = -2 * t379;
t163 = -0.2e1 * t377;
t164 = 0.2e1 * t378;
t166 = mrSges(3,1) * t390;
t221 = pkin(1) ^ 2;
t234 = 0.2e1 * mrSges(3,3) * pkin(5) + (t200 * t221) + Ifges(2,1) + Ifges(3,2) + Ifges(1,3);
t220 = pkin(2) ^ 2;
t250 = t220 * m(3) - Ifges(2,1) + Ifges(2,2) + t179;
t362 = mrSges(3,2) * t171;
t289 = pkin(1) * t362;
t358 = Ifges(3,4) * t184;
t295 = 0.4e1 * t358;
t299 = pkin(1) * t385;
t364 = (2 * Ifges(2,4)) - 0.2e1 * Ifges(3,4);
t31 = (-0.2e1 * t315 + (t164 + t295) * t190 - 0.2e1 * t288 + t250) * t156 + (t166 * t190 + (t337 - t361) * t390 + (0.4e1 * t354 + t163 * t190 + 0.2e1 * (t179 * t190 - t378) * t184 + t364) * t171) * t172 + t315 + 0.2e1 * (-t289 - t358) * t190 + (t136 + mrSges(2,2)) * t299 + (t176 ^ 2 * m(3)) + t234 + ((t384 + t401) * qJ(2,3));
t312 = t179 * t184;
t43 = t55 * t312;
t58 = (-t103 * t190 - t106 * t184) * t172 + t171 * (t103 * t184 - t106 * t190);
t61 = -t149 + t246;
t210 = 3 * pkin(7);
t214 = 2 * qJ(3,3);
t217 = pkin(3) ^ 2;
t254 = pkin(3) * t306 / 0.2e1;
t266 = -0.3e1 * t217 - 0.2e1 * t220 - (4 * t221);
t267 = -0.4e1 * pkin(3) * t129;
t304 = -0.2e1 * pkin(2) * pkin(3);
t7 = (t110 + t373) * t389 * t394 * t336 + (t389 * t376 * t400 + (-(-t217 * cos((3 * t160)) + (-0.4e1 * t150 ^ 2 + t266) * t126 + t267 + (cos((t214 + t210)) + cos((t214 + pkin(7)))) * t304 + (-cos((-pkin(7) + qJ(3,3))) - cos((t210 + qJ(3,3)))) * t220) * t55 / 0.4e1 - t114 * t254 * t400 + (-t410 - t413) * (-t53 - (-t343 / 0.2e1 - t340 / 0.2e1 - t277 / 0.4e1) * t142) - (t126 * t390 + t407) * t28) * t114) * t142 * t55;
t82 = -g(1) * t130 + g(2) * t133;
t85 = g(1) * t133 + g(2) * t130;
t278 = m(3) * pkin(5) + t180;
t91 = (qJ(2,3) * t200) + t278;
t1 = t31 * t10 - t61 * t7 + t58 * t253 - 0.2e1 * ((((-qJ(2,3) / 0.2e1 + t213) * mrSges(3,2) + t201) * t273 - t298) * t171 - t43) * t190 * t273 - 0.2e1 * t55 * (Ifges(3,4) * t273 - t91 * t28) + (-t235 * t85 + t240 * t82) * t191 - (t82 * t235 + t240 * t85) * t185 - 0.4e1 * (((t147 + t312) * t190 + t184 * t280 + t395) * t156 - t354) * t55 * t273 + ((-t258 + (((-qJ(2,3) / 0.4e1 + t212) * mrSges(3,1) + t202) * t273 - ((t148 + t358) * t385 + t281) * t55) * t190 - t284 / 0.2e1 - t184 * (-t103 * t273 - 0.2e1 * t298) / 0.4e1) * t261 - (-t106 * t273 - t55 * t300) * t184 * t269) * t114;
t370 = t142 * t1;
t262 = 0.2e1 * t306;
t4 = t61 * t10 + (-t91 * t55 + (-t249 * t171 - t303 * t172) * t114 * t262) * t55 + (-t185 * t85 + t191 * t82 - t7) * t200;
t369 = t142 * t4;
t117 = t116 * t388;
t276 = t65 * t331;
t68 = t110 * t187 + t151 * t193;
t71 = t110 * t193 - t151 * t187;
t50 = t131 * t71 + t134 * t68 + t77 * t372;
t339 = t50 * t198;
t47 = -t131 * t68 + t134 * t71 + t74 * t372;
t342 = t47 * t199;
t229 = t276 + 0.2e1 * t339 + 0.2e1 * t342;
t333 = t116 * t143;
t264 = -t333 / 0.2e1;
t29 = (t342 + t339 + t276 / 0.2e1) * t143;
t11 = (-t117 * t394 - (-t29 + t116 * t56 * t408 + (t54 * t116 + t229 * t264) * t127) * t56) * t143;
t236 = m(3) * t177 + t285 + t402;
t360 = mrSges(3,2) * t186;
t248 = mrSges(3,1) * t192 - t360;
t137 = t186 * mrSges(3,1);
t302 = mrSges(3,2) * t192 + t137;
t245 = t171 * (mrSges(2,2) + t302) + t172 * (-t248 - t337);
t239 = -t109 + t245;
t252 = t117 * t121 * t317;
t314 = t179 * t168;
t350 = t171 * t56;
t257 = t314 * t350;
t271 = t116 * t306;
t287 = pkin(2) * t360;
t283 = (t287 - t179) * t350;
t297 = mrSges(3,1) * t54;
t311 = t179 * t186;
t357 = Ifges(3,4) * t186;
t294 = 0.4e1 * t357;
t32 = (-0.2e1 * t314 + (t164 + t294) * t192 - 0.2e1 * t287 + t250) * t156 + (t166 * t192 + (t337 - t360) * t390 + (0.4e1 * t353 + t163 * t192 + 0.2e1 * (t179 * t192 - t378) * t186 + t364) * t171) * t172 + t314 + 0.2e1 * (-t289 - t357) * t192 + (t137 + mrSges(2,2)) * t299 + (t177 ^ 2 * m(3)) + t234 + ((t384 + t402) * qJ(2,2));
t44 = t56 * t311;
t59 = (-t104 * t192 - t107 * t186) * t172 + t171 * (t104 * t186 - t107 * t192);
t62 = -t149 + t245;
t215 = 2 * qJ(3,2);
t8 = (t110 + t372) * t388 * t394 * t333 + (t388 * t376 * t399 + (-(-t217 * cos((3 * t161)) + (-0.4e1 * t151 ^ 2 + t266) * t127 + t267 + (cos((t210 + t215)) + cos((t215 + pkin(7)))) * t304 + (-cos((qJ(3,2) + t210)) - cos((-pkin(7) + qJ(3,2)))) * t220) * t56 / 0.4e1 - t116 * t254 * t399 + (-t411 - t414) * (-t54 - (-t342 / 0.2e1 - t339 / 0.2e1 - t276 / 0.4e1) * t143) - (t127 * t390 + t408) * t29) * t116) * t143 * t56;
t83 = -g(1) * t131 + g(2) * t134;
t86 = g(1) * t134 + g(2) * t131;
t92 = (qJ(2,2) * t200) + t278;
t2 = t32 * t11 - t62 * t8 + t59 * t252 - 0.2e1 * ((((-qJ(2,2) / 0.2e1 + t213) * mrSges(3,2) + t201) * t271 - t297) * t171 - t44) * t192 * t271 - 0.2e1 * t56 * (Ifges(3,4) * t271 - t92 * t29) + (-t236 * t86 + t239 * t83) * t193 - (t83 * t236 + t239 * t86) * t187 - 0.4e1 * (((t147 + t311) * t192 + t186 * t280 + t396) * t156 - t353) * t56 * t271 + ((-t257 + (((-qJ(2,2) / 0.4e1 + t212) * mrSges(3,1) + t202) * t271 - ((t148 + t357) * t385 + t281) * t56) * t192 - t283 / 0.2e1 - t186 * (-t104 * t271 - 0.2e1 * t297) / 0.4e1) * t261 - (-t107 * t271 - t56 * t300) * t186 * t269) * t116;
t368 = t143 * t2;
t5 = t62 * t11 + (-t92 * t56 + (-t248 * t171 - t302 * t172) * t116 * t262) * t56 + (-t187 * t86 + t193 * t83 - t8) * t200;
t367 = t143 * t5;
t119 = t118 * t387;
t330 = t118 * t144;
t263 = -t330 / 0.2e1;
t30 = (t341 + t338 + t275 / 0.2e1) * t144;
t12 = (-t119 * t394 - (-t30 + t118 * t57 * t409 + (t52 * t118 + t228 * t263) * t128) * t57) * t144;
t237 = m(3) * t178 + t285 + t403;
t359 = mrSges(3,2) * t188;
t247 = mrSges(3,1) * t194 - t359;
t138 = t188 * mrSges(3,1);
t301 = mrSges(3,2) * t194 + t138;
t244 = t171 * (mrSges(2,2) + t301) + t172 * (-t247 - t337);
t238 = -t109 + t244;
t251 = t119 * t122 * t317;
t313 = t179 * t169;
t349 = t171 * t57;
t256 = t313 * t349;
t268 = t118 * t306;
t286 = pkin(2) * t359;
t282 = (t286 - t179) * t349;
t296 = mrSges(3,1) * t52;
t310 = t179 * t188;
t316 = t171 * t188;
t356 = Ifges(3,4) * t188;
t293 = 0.4e1 * t356;
t33 = (-0.2e1 * t313 + (t164 + t293) * t194 - 0.2e1 * t286 + t250) * t156 + (t166 * t194 + (t337 - t359) * t390 + (0.4e1 * t352 + t163 * t194 + 0.2e1 * (t179 * t194 - t378) * t188 + t364) * t171) * t172 + t313 + 0.2e1 * (-t289 - t356) * t194 + (t138 + mrSges(2,2)) * t299 + (t178 ^ 2 * m(3)) + t234 + ((t384 + t403) * qJ(2,1));
t45 = t57 * t310;
t60 = (-t105 * t194 - t108 * t188) * t172 + t171 * (t105 * t188 - t108 * t194);
t63 = -t149 + t244;
t84 = -g(1) * t132 + g(2) * t135;
t87 = g(1) * t135 + g(2) * t132;
t216 = 2 * qJ(3,1);
t9 = (t110 + t371) * t387 * t394 * t330 + (t387 * t376 * t398 + (-(-t217 * cos((3 * t162)) + (-0.4e1 * t152 ^ 2 + t266) * t128 + t267 + (cos((t210 + t216)) + cos((pkin(7) + t216))) * t304 + (-cos((-pkin(7) + qJ(3,1))) - cos((t210 + qJ(3,1)))) * t220) * t57 / 0.4e1 - t118 * t254 * t398 + (-t412 - t415) * (-t52 - (-t341 / 0.2e1 - t338 / 0.2e1 - t275 / 0.4e1) * t144) - (t128 * t390 + t409) * t30) * t118) * t144 * t57;
t93 = (qJ(2,1) * t200) + t278;
t3 = t33 * t12 - t63 * t9 + t60 * t251 + t118 * (-t256 + (((-qJ(2,1) / 0.4e1 + t212) * mrSges(3,1) + t202) * t268 - ((t148 + t356) * t385 + t281) * t57) * t194 - t282 / 0.2e1 - t188 * (-t105 * t268 - 0.2e1 * t296) / 0.4e1) * t261 - 0.2e1 * ((((-qJ(2,1) / 0.2e1 + t213) * mrSges(3,2) + t201) * t268 - t296) * t171 - t45) * t194 * t268 - (-t108 * t268 - t57 * t300) * t268 * t316 - 0.2e1 * t57 * (Ifges(3,4) * t268 - t93 * t30) + (-t237 * t87 + t238 * t84) * t195 - (t84 * t237 + t238 * t87) * t189 - 0.4e1 * (((t147 + t310) * t194 + t188 * t280 + t397) * t156 - t352) * t57 * t268;
t366 = t144 * t3;
t6 = t63 * t12 + (-t93 * t57 + (-t247 * t171 - t301 * t172) * t118 * t262) * t57 + (-t189 * t87 + t195 * t84 - t9) * t200;
t365 = t144 * t6;
t363 = mrSges(3,1) * t171;
t355 = Ifges(3,3) * t218;
t26 = t230 * t142 - t53;
t348 = t184 * t26;
t27 = t229 * t143 - t54;
t347 = t186 * t27;
t346 = t218 * t58;
t345 = t218 * t59;
t344 = t218 * t60;
t181 = xDDP(3);
t335 = t114 * t181;
t332 = t116 * t181;
t329 = t118 * t181;
t327 = t120 * t142;
t326 = t121 * t143;
t325 = t122 * t144;
t182 = xDDP(2);
t324 = t142 * t182;
t183 = xDDP(1);
t323 = t142 * t183;
t322 = t143 * t182;
t321 = t143 * t183;
t320 = t144 * t182;
t319 = t144 * t183;
t309 = t218 * t114;
t308 = t218 * t116;
t307 = t218 * t118;
t279 = -t378 / 0.4e1;
t274 = t142 * t335;
t272 = t143 * t332;
t270 = t144 * t329;
t243 = t185 * t82 + t191 * t85;
t242 = t187 * t83 + t193 * t86;
t241 = t189 * t84 + t195 * t87;
t203 = -Ifges(3,4) / 0.2e1;
t196 = mrSges(3,2) * g(3);
t146 = -t377 / 0.4e1;
t145 = Ifges(3,1) / 0.2e1 - Ifges(3,2) / 0.2e1;
t42 = (t122 * t63 + t66 * t375) * t330;
t41 = (t121 * t62 + t65 * t375) * t333;
t40 = (t120 * t61 + t64 * t375) * t336;
t39 = (t200 * t51 + t63 * t78) * t144;
t38 = (t200 * t50 + t62 * t77) * t143;
t37 = (t200 * t49 + t61 * t76) * t142;
t36 = (t200 * t48 + t63 * t75) * t144;
t35 = (t200 * t47 + t62 * t74) * t143;
t34 = (t200 * t46 + t61 * t73) * t142;
t21 = (t33 * t78 + t51 * t63) * t144;
t20 = (t32 * t77 + t50 * t62) * t143;
t19 = (t31 * t76 + t49 * t61) * t142;
t18 = (t33 * t75 + t48 * t63) * t144;
t17 = (t32 * t74 + t47 * t62) * t143;
t16 = (t31 * t73 + t46 * t61) * t142;
t15 = (t344 - (t122 * t33 + t63 * t381) * t144) * t118;
t14 = (t345 - (t121 * t32 + t62 * t382) * t143) * t116;
t13 = (t346 - (t120 * t31 + t61 * t383) * t142) * t114;
t22 = [-(-t18 * t75 - t36 * t48) * t319 - (-t18 * t78 - t36 * t51) * t320 - (-t18 * t122 + t75 * t344 - t36 * t381) * t270 - t75 * t366 - t48 * t365 - (-t17 * t74 - t35 * t47) * t321 - (-t17 * t77 - t35 * t50) * t322 - (-t17 * t121 + t74 * t345 - t35 * t382) * t272 - t74 * t368 - t47 * t367 - (-t16 * t73 - t34 * t46) * t323 - (-t16 * t76 - t34 * t49) * t324 - (-t16 * t120 - t34 * t383 + t73 * t346) * t274 - t73 * t370 - t46 * t369 + (-g(1) + t183) * m(4); -(-t21 * t75 - t39 * t48) * t319 - (-t21 * t78 - t39 * t51) * t320 - (-t21 * t122 + t78 * t344 - t39 * t381) * t270 - t78 * t366 - t51 * t365 - (-t20 * t74 - t38 * t47) * t321 - (-t20 * t77 - t38 * t50) * t322 - (-t20 * t121 + t77 * t345 - t38 * t382) * t272 - t77 * t368 - t50 * t367 - (-t19 * t73 - t37 * t46) * t323 - (-t19 * t76 - t37 * t49) * t324 - (-t19 * t120 + t76 * t346 - t37 * t383) * t274 - t76 * t370 - t49 * t369 + (-g(2) + t182) * m(4); -(t15 * t75 - t42 * t48) * t319 - (t15 * t78 - t42 * t51) * t320 + (-t15 * t325 + t42 * t144 * t381 + (-t60 * t325 + t355) * t307) * t329 - t118 * t3 * t325 + t66 * t6 * t263 + (t60 * t12 + Ifges(3,3) * t251 + (t241 * mrSges(3,2) - t380) * t128 + t122 * (t241 * mrSges(3,1) + t196) - ((0.2e1 * t256 + (t374 - (t293 + t378) * t349) * t194 + t282 + t25 * t138) * t172 + (t25 * t363 - t45) * t194 - t316 * t374 - ((t352 + (t145 * t188 + t146) * t194 + t188 * t279 + t203) * t386 + t397) * t57) * t57) * t307 - (t14 * t74 - t41 * t47) * t321 - (t14 * t77 - t41 * t50) * t322 + (-t14 * t326 + t41 * t143 * t382 + (-t59 * t326 + t355) * t308) * t332 - t116 * t2 * t326 + t65 * t5 * t264 + (t59 * t11 + Ifges(3,3) * t252 + (t242 * mrSges(3,2) - t380) * t127 + t121 * (t242 * mrSges(3,1) + t196) - ((0.2e1 * t257 + (mrSges(3,2) * t27 - (t294 + t378) * t350) * t192 + t283 + mrSges(3,1) * t347) * t172 + (t27 * t363 - t44) * t192 - t347 * t362 - ((t353 + (t145 * t186 + t146) * t192 + t186 * t279 + t203) * t386 + t396) * t56) * t56) * t308 - (t13 * t73 - t40 * t46) * t323 - (t13 * t76 - t40 * t49) * t324 + (-t13 * t327 + t40 * t142 * t383 + (-t58 * t327 + t355) * t309) * t335 - t114 * t1 * t327 + t64 * t4 * t265 + (t58 * t10 + Ifges(3,3) * t253 + (t243 * mrSges(3,2) - t380) * t126 + t120 * (t243 * mrSges(3,1) + t196) - ((0.2e1 * t258 + (mrSges(3,2) * t26 - (t295 + t378) * t351) * t190 + t284 + mrSges(3,1) * t348) * t172 + (t26 * t363 - t43) * t190 - t348 * t362 - ((t354 + (t145 * t184 + t146) * t190 + t184 * t279 + t203) * t386 + t395) * t55) * t55) * t309 + (-g(3) + t181) * m(4);];
tauX  = t22;
