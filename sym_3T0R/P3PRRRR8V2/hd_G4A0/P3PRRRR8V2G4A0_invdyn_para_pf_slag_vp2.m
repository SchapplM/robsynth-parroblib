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
% Datum: 2020-08-06 18:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRRR8V2G4A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G4A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V2G4A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR8V2G4A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G4A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V2G4A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G4A0_invdyn_para_pf_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G4A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V2G4A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR8V2G4A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G4A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G4A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:14:41
% EndTime: 2020-08-06 18:14:55
% DurationCPUTime: 13.02s
% Computational Cost: add. (62559->577), mult. (127572->1075), div. (5112->10), fcn. (147594->34), ass. (0->421)
t273 = cos(qJ(2,1));
t279 = pkin(7) + pkin(6);
t210 = t273 * t279;
t267 = sin(qJ(2,1));
t189 = pkin(2) * t267 - t210;
t246 = sin(pkin(4));
t248 = cos(pkin(4));
t266 = sin(qJ(3,1));
t356 = t248 * t266;
t153 = pkin(3) * t356 + t189 * t246;
t272 = cos(qJ(3,1));
t367 = t246 * t267;
t244 = t272 ^ 2;
t456 = pkin(3) * t244;
t120 = 0.1e1 / (pkin(2) * t356 + t153 * t272 + t367 * t456);
t271 = cos(qJ(2,2));
t209 = t271 * t279;
t265 = sin(qJ(2,2));
t188 = pkin(2) * t265 - t209;
t264 = sin(qJ(3,2));
t358 = t248 * t264;
t152 = pkin(3) * t358 + t188 * t246;
t270 = cos(qJ(3,2));
t369 = t246 * t265;
t243 = t270 ^ 2;
t457 = pkin(3) * t243;
t119 = 0.1e1 / (pkin(2) * t358 + t152 * t270 + t369 * t457);
t269 = cos(qJ(2,3));
t208 = t269 * t279;
t263 = sin(qJ(2,3));
t187 = pkin(2) * t263 - t208;
t262 = sin(qJ(3,3));
t360 = t248 * t262;
t151 = pkin(3) * t360 + t187 * t246;
t268 = cos(qJ(3,3));
t371 = t246 * t263;
t242 = t268 ^ 2;
t458 = pkin(3) * t242;
t118 = 0.1e1 / (pkin(2) * t360 + t151 * t268 + t371 * t458);
t249 = legFrame(3,3);
t213 = sin(t249);
t219 = cos(t249);
t245 = sin(pkin(8));
t247 = cos(pkin(8));
t157 = -t213 * t245 + t219 * t247;
t160 = t213 * t247 + t219 * t245;
t205 = pkin(3) * t268 + pkin(2);
t178 = t205 * t263 - t208;
t349 = t263 * t279;
t387 = (t205 * t269 + t349) * t248;
t100 = -t157 * t387 + t160 * t178;
t276 = xDP(3);
t277 = xDP(2);
t256 = legFrame(3,2);
t233 = cos(t256);
t278 = xDP(1);
t377 = t233 * t278;
t184 = t205 * t360;
t366 = t246 * t268;
t136 = 0.1e1 / (t178 * t366 + t184);
t283 = 0.1e1 / pkin(3);
t402 = t136 * t283;
t252 = legFrame(3,1);
t222 = cos(t252);
t216 = sin(t252);
t230 = sin(t256);
t384 = t216 * t230;
t115 = t157 * t384 + t160 * t222;
t82 = -t115 * t387 + (-t157 * t222 + t160 * t384) * t178;
t381 = t222 * t230;
t112 = t157 * t381 - t160 * t216;
t85 = t112 * t387 - (t157 * t216 + t160 * t381) * t178;
t64 = (t100 * t377 + t276 * t85 + t277 * t82) * t402;
t480 = -0.2e1 * t64;
t250 = legFrame(2,3);
t214 = sin(t250);
t220 = cos(t250);
t158 = -t214 * t245 + t220 * t247;
t161 = t214 * t247 + t220 * t245;
t206 = pkin(3) * t270 + pkin(2);
t179 = t206 * t265 - t209;
t347 = t265 * t279;
t386 = (t206 * t271 + t347) * t248;
t101 = -t158 * t386 + t161 * t179;
t257 = legFrame(2,2);
t234 = cos(t257);
t375 = t234 * t278;
t185 = t206 * t358;
t364 = t246 * t270;
t137 = 0.1e1 / (t179 * t364 + t185);
t401 = t137 * t283;
t253 = legFrame(2,1);
t223 = cos(t253);
t217 = sin(t253);
t231 = sin(t257);
t383 = t217 * t231;
t116 = t158 * t383 + t161 * t223;
t83 = -t116 * t386 + (-t158 * t223 + t161 * t383) * t179;
t380 = t223 * t231;
t113 = t158 * t380 - t161 * t217;
t86 = t113 * t386 - (t158 * t217 + t161 * t380) * t179;
t65 = (t101 * t375 + t276 * t86 + t277 * t83) * t401;
t479 = -0.2e1 * t65;
t251 = legFrame(1,3);
t215 = sin(t251);
t221 = cos(t251);
t159 = -t215 * t245 + t221 * t247;
t162 = t215 * t247 + t221 * t245;
t207 = pkin(3) * t272 + pkin(2);
t180 = t207 * t267 - t210;
t345 = t267 * t279;
t385 = (t207 * t273 + t345) * t248;
t102 = -t159 * t385 + t162 * t180;
t258 = legFrame(1,2);
t235 = cos(t258);
t373 = t235 * t278;
t186 = t207 * t356;
t362 = t246 * t272;
t138 = 0.1e1 / (t180 * t362 + t186);
t400 = t138 * t283;
t254 = legFrame(1,1);
t224 = cos(t254);
t218 = sin(t254);
t232 = sin(t258);
t382 = t218 * t232;
t117 = t159 * t382 + t162 * t224;
t84 = -t117 * t385 + (-t159 * t224 + t162 * t382) * t180;
t379 = t224 * t232;
t114 = t159 * t379 - t162 * t218;
t87 = t114 * t385 - (t159 * t218 + t162 * t379) * t180;
t66 = (t102 * t373 + t276 * t87 + t277 * t84) * t400;
t478 = -0.2e1 * t66;
t312 = t268 * mrSges(3,1) - mrSges(3,2) * t262;
t311 = t270 * mrSges(3,1) - mrSges(3,2) * t264;
t310 = t272 * mrSges(3,1) - mrSges(3,2) * t266;
t147 = g(1) * t232 + (-g(2) * t218 + g(3) * t224) * t235;
t453 = g(1) * t235;
t110 = -t215 * t453 + (-t215 * t382 + t221 * t224) * g(2) + (t215 * t379 + t218 * t221) * g(3);
t111 = t221 * t453 + (t215 * t224 + t221 * t382) * g(2) + (t215 * t218 - t221 * t379) * g(3);
t296 = t110 * t247 - t111 * t245;
t474 = -t147 * t248 + t296 * t246;
t146 = g(1) * t231 + (-g(2) * t217 + g(3) * t223) * t234;
t454 = g(1) * t234;
t108 = -t214 * t454 + (-t214 * t383 + t220 * t223) * g(2) + (t214 * t380 + t217 * t220) * g(3);
t109 = t220 * t454 + (t214 * t223 + t220 * t383) * g(2) + (t214 * t217 - t220 * t380) * g(3);
t298 = t108 * t247 - t109 * t245;
t473 = -t146 * t248 + t298 * t246;
t145 = g(1) * t230 + (-g(2) * t216 + g(3) * t222) * t233;
t455 = g(1) * t233;
t106 = -t213 * t455 + (-t213 * t384 + t219 * t222) * g(2) + (t213 * t381 + t216 * t219) * g(3);
t107 = t219 * t455 + (t213 * t222 + t219 * t384) * g(2) + (t213 * t216 - t219 * t381) * g(3);
t300 = t106 * t247 - t107 * t245;
t472 = -t145 * t248 + t300 * t246;
t471 = t147 * t246 + t296 * t248;
t470 = t146 * t246 + t298 * t248;
t469 = t145 * t246 + t300 * t248;
t255 = Ifges(3,1) - Ifges(3,2);
t274 = mrSges(3,2) * pkin(2);
t465 = -2 * Ifges(3,4);
t468 = Ifges(3,4) + t244 * t465 + (-t255 * t266 + t274) * t272;
t467 = Ifges(3,4) + t243 * t465 + (-t255 * t264 + t274) * t270;
t466 = Ifges(3,4) + t242 * t465 + (-t255 * t262 + t274) * t268;
t275 = mrSges(3,1) * pkin(2);
t464 = pkin(3) * t64;
t463 = pkin(3) * t65;
t462 = pkin(3) * t66;
t225 = pkin(6) * mrSges(3,2) - Ifges(3,6);
t461 = -t225 / 0.2e1;
t226 = mrSges(3,1) * pkin(6) - Ifges(3,5);
t460 = t226 / 0.2e1;
t459 = pkin(2) * t262;
t284 = pkin(2) ^ 2;
t211 = t279 ^ 2 + t284;
t282 = pkin(3) ^ 2;
t336 = t262 * t464;
t344 = 0.2e1 * pkin(2) * pkin(3);
t350 = t263 * t268;
t154 = pkin(3) * t350 + t187;
t133 = 0.1e1 / (t154 * t366 + t184);
t359 = t248 * t263;
t163 = t245 * t359 - t247 * t269;
t166 = t245 * t269 + t247 * t359;
t121 = -t163 * t213 + t166 * t219;
t294 = t163 * t219 + t166 * t213;
t76 = (t121 * t381 - t216 * t294) * t262 + t112 * t366;
t79 = (-t121 * t384 - t294 * t222) * t262 - t115 * t366;
t97 = t157 * t366 + t262 * (t157 * t359 + t160 * t269);
t58 = -t118 * t97 * t377 + (t276 * t76 + t277 * t79) * t133;
t452 = (-t279 * t336 + (t242 * t282 + t268 * t344 + t211) * t58) * t58;
t335 = t264 * t463;
t348 = t265 * t270;
t155 = pkin(3) * t348 + t188;
t134 = 0.1e1 / (t155 * t364 + t185);
t357 = t248 * t265;
t164 = t245 * t357 - t247 * t271;
t167 = t245 * t271 + t247 * t357;
t122 = -t164 * t214 + t167 * t220;
t293 = t164 * t220 + t167 * t214;
t77 = (t122 * t380 - t217 * t293) * t264 + t113 * t364;
t80 = (-t122 * t383 - t293 * t223) * t264 - t116 * t364;
t98 = t158 * t364 + t264 * (t158 * t357 + t161 * t271);
t59 = -t119 * t98 * t375 + (t276 * t77 + t277 * t80) * t134;
t451 = (-t279 * t335 + (t243 * t282 + t270 * t344 + t211) * t59) * t59;
t334 = t266 * t462;
t346 = t267 * t272;
t156 = pkin(3) * t346 + t189;
t135 = 0.1e1 / (t156 * t362 + t186);
t355 = t248 * t267;
t165 = t245 * t355 - t247 * t273;
t168 = t245 * t273 + t247 * t355;
t123 = -t165 * t215 + t168 * t221;
t292 = t165 * t221 + t168 * t215;
t78 = (t123 * t379 - t218 * t292) * t266 + t114 * t362;
t81 = (-t123 * t382 - t292 * t224) * t266 - t117 * t362;
t99 = t159 * t362 + t266 * (t159 * t355 + t162 * t273);
t60 = -t120 * t99 * t373 + (t276 * t78 + t277 * t81) * t135;
t450 = (-t279 * t334 + (t244 * t282 + t272 * t344 + t211) * t60) * t60;
t449 = t264 * pkin(2);
t448 = t266 * pkin(2);
t447 = mrSges(3,1) * t262;
t446 = mrSges(3,1) * t264;
t445 = mrSges(3,1) * t266;
t103 = t294 * t230 + t233 * t371;
t130 = -t160 * t246 * t230 + t233 * t248;
t390 = t157 * t246;
t190 = pkin(2) * t269 + t349;
t372 = t246 * t262;
t291 = pkin(3) * t372 - t187 * t248;
t124 = -t190 * t247 - t291 * t245;
t127 = t190 * t245 - t291 * t247;
t88 = t151 * t233 + (t124 * t219 + t127 * t213) * t230;
t94 = -t124 * t213 + t127 * t219;
t67 = -(t103 * t216 - t121 * t222) * t458 + (-t216 * t88 + t222 * t94) * t268 - (t130 * t216 + t222 * t390) * t459;
t441 = t118 * t67;
t68 = (t103 * t222 + t121 * t216) * t458 + (t216 * t94 + t222 * t88) * t268 + (t130 * t222 - t216 * t390) * t459;
t440 = t118 * t68;
t104 = t293 * t231 + t234 * t369;
t131 = -t161 * t246 * t231 + t234 * t248;
t389 = t158 * t246;
t191 = pkin(2) * t271 + t347;
t370 = t246 * t264;
t290 = pkin(3) * t370 - t188 * t248;
t125 = -t191 * t247 - t290 * t245;
t128 = t191 * t245 - t290 * t247;
t89 = t152 * t234 + (t125 * t220 + t128 * t214) * t231;
t95 = -t125 * t214 + t128 * t220;
t69 = -(t104 * t217 - t122 * t223) * t457 + (-t217 * t89 + t223 * t95) * t270 - (t131 * t217 + t223 * t389) * t449;
t439 = t119 * t69;
t70 = (t104 * t223 + t122 * t217) * t457 + (t217 * t95 + t223 * t89) * t270 + (t131 * t223 - t217 * t389) * t449;
t438 = t119 * t70;
t105 = t292 * t232 + t235 * t367;
t132 = -t162 * t246 * t232 + t235 * t248;
t388 = t159 * t246;
t192 = pkin(2) * t273 + t345;
t368 = t246 * t266;
t289 = pkin(3) * t368 - t189 * t248;
t126 = -t192 * t247 - t289 * t245;
t129 = t192 * t245 - t289 * t247;
t90 = t153 * t235 + (t126 * t221 + t129 * t215) * t232;
t96 = -t126 * t215 + t129 * t221;
t71 = -(t105 * t218 - t123 * t224) * t456 + (-t218 * t90 + t224 * t96) * t272 - (t132 * t218 + t224 * t388) * t448;
t437 = t120 * t71;
t72 = (t105 * t224 + t123 * t218) * t456 + (t218 * t96 + t224 * t90) * t272 + (t132 * t224 - t218 * t388) * t448;
t436 = t120 * t72;
t435 = t133 * t76;
t434 = t133 * t79;
t433 = t134 * t77;
t432 = t134 * t80;
t431 = t135 * t78;
t430 = t135 * t81;
t429 = t136 * t82;
t428 = t136 * t85;
t427 = t137 * t83;
t426 = t137 * t86;
t425 = t138 * t84;
t424 = t138 * t87;
t423 = t233 * t97;
t422 = t234 * t98;
t421 = t235 * t99;
t420 = t269 * t58;
t419 = t271 * t59;
t418 = t273 * t60;
t417 = t279 * t58;
t416 = t279 * t59;
t415 = t279 * t60;
t236 = m(3) * pkin(2) + mrSges(2,1);
t306 = mrSges(3,2) * t268 + t447;
t142 = t312 * t248 - t306 * t371;
t414 = t118 * t142;
t240 = m(1) + m(2) + m(3);
t413 = t118 * t240;
t305 = mrSges(3,2) * t270 + t446;
t143 = t311 * t248 - t305 * t369;
t412 = t119 * t143;
t411 = t119 * t240;
t304 = mrSges(3,2) * t272 + t445;
t144 = t310 * t248 - t304 * t367;
t410 = t120 * t144;
t409 = t120 * t240;
t288 = Ifges(3,1) + Ifges(2,3) + (pkin(6) ^ 2 + t284) * m(3) + 0.2e1 * pkin(6) * mrSges(3,3);
t340 = -0.2e1 * t274;
t139 = -t255 * t242 + 0.2e1 * (Ifges(3,4) * t262 + t275) * t268 + t262 * t340 + t288;
t408 = t133 * t139;
t169 = -t225 * t268 - t262 * t226;
t407 = t133 * t169;
t140 = -t255 * t243 + 0.2e1 * (Ifges(3,4) * t264 + t275) * t270 + t264 * t340 + t288;
t406 = t134 * t140;
t170 = -t225 * t270 - t264 * t226;
t405 = t134 * t170;
t141 = -t255 * t244 + 0.2e1 * (Ifges(3,4) * t266 + t275) * t272 + t266 * t340 + t288;
t404 = t135 * t141;
t171 = -t225 * t272 - t266 * t226;
t403 = t135 * t171;
t181 = t312 + t236;
t212 = m(3) * pkin(6) - mrSges(2,2) + mrSges(3,3);
t148 = t181 * t269 + t212 * t263;
t393 = t148 * t246;
t182 = t311 + t236;
t149 = t182 * t271 + t212 * t265;
t392 = t149 * t246;
t183 = t310 + t236;
t150 = t183 * t273 + t212 * t267;
t391 = t150 * t246;
t378 = t233 * t246;
t376 = t234 * t246;
t374 = t235 * t246;
t365 = t246 * t269;
t363 = t246 * t271;
t361 = t246 * t273;
t354 = t248 * t283;
t333 = Ifges(3,3) * t402;
t332 = Ifges(3,3) * t401;
t331 = Ifges(3,3) * t400;
t330 = t262 * t417;
t329 = t264 * t416;
t328 = t266 * t415;
t327 = t100 * t136 * t233;
t326 = t101 * t137 * t234;
t325 = t102 * t138 * t235;
t324 = t118 * t393;
t323 = t119 * t392;
t322 = t120 * t391;
t321 = t133 * t393;
t320 = t134 * t392;
t319 = t135 * t391;
t318 = t142 * t402;
t317 = t169 * t402;
t316 = t143 * t401;
t315 = t170 * t401;
t314 = t144 * t400;
t313 = t171 * t400;
t309 = t283 * t327;
t308 = t283 * t326;
t307 = t283 * t325;
t299 = t106 * t245 + t107 * t247;
t297 = t108 * t245 + t109 * t247;
t295 = t110 * t245 + t111 * t247;
t287 = t469 * t263 + t299 * t269;
t286 = t470 * t265 + t297 * t271;
t285 = t471 * t267 + t295 * t273;
t261 = xDDP(1);
t260 = xDDP(2);
t259 = xDDP(3);
t75 = -((-t159 * t273 + t162 * t355) * t235 - t232 * t367) * t456 + ((t159 * t192 + t289 * t162) * t235 + t153 * t232) * t272 + (t162 * t374 + t232 * t248) * t448;
t74 = -((-t158 * t271 + t161 * t357) * t234 - t231 * t369) * t457 + ((t158 * t191 + t290 * t161) * t234 + t152 * t231) * t270 + (t161 * t376 + t231 * t248) * t449;
t73 = -((-t157 * t269 + t160 * t359) * t233 - t230 * t371) * t458 + ((t157 * t190 + t291 * t160) * t233 + t151 * t230) * t268 + (t160 * t378 + t230 * t248) * t459;
t63 = t66 ^ 2;
t62 = t65 ^ 2;
t61 = t64 ^ 2;
t57 = t60 ^ 2;
t56 = t59 ^ 2;
t55 = t58 ^ 2;
t54 = Ifges(3,3) * t307 + (t144 * t75 - t171 * t421) * t120;
t53 = Ifges(3,3) * t308 + (t143 * t74 - t170 * t422) * t119;
t52 = Ifges(3,3) * t309 + (t142 * t73 - t169 * t423) * t118;
t51 = t144 * t307 + (-t150 * t99 * t374 + t240 * t75) * t120;
t50 = t143 * t308 + (-t149 * t98 * t376 + t240 * t74) * t119;
t49 = t142 * t309 + (-t148 * t97 * t378 + t240 * t73) * t118;
t48 = t171 * t307 + (-t141 * t421 + t75 * t391) * t120;
t47 = t170 * t308 + (-t140 * t422 + t74 * t392) * t119;
t46 = t169 * t309 + (-t139 * t423 + t73 * t393) * t118;
t45 = t87 * t331 + t78 * t403 + t72 * t410;
t44 = t84 * t331 + t81 * t403 + t71 * t410;
t43 = t86 * t332 + t77 * t405 + t70 * t412;
t42 = t83 * t332 + t80 * t405 + t69 * t412;
t41 = t85 * t333 + t76 * t407 + t68 * t414;
t40 = t82 * t333 + t79 * t407 + t67 * t414;
t39 = t87 * t314 + t78 * t319 + t72 * t409;
t38 = t84 * t314 + t81 * t319 + t71 * t409;
t37 = t86 * t316 + t77 * t320 + t70 * t411;
t36 = t83 * t316 + t80 * t320 + t69 * t411;
t35 = t85 * t318 + t76 * t321 + t68 * t413;
t34 = t82 * t318 + t79 * t321 + t67 * t413;
t33 = t87 * t313 + t72 * t322 + t78 * t404;
t32 = t84 * t313 + t71 * t322 + t81 * t404;
t31 = t86 * t315 + t70 * t323 + t77 * t406;
t30 = t83 * t315 + t69 * t323 + t80 * t406;
t29 = t85 * t317 + t68 * t324 + t76 * t408;
t28 = t82 * t317 + t67 * t324 + t79 * t408;
t24 = t328 - t462;
t23 = t329 - t463;
t22 = t330 - t464;
t18 = (-t272 * t450 - (pkin(2) * t66 - t24 * t272) * t462) * t120;
t17 = (-t270 * t451 - (pkin(2) * t65 - t23 * t270) * t463) * t119;
t16 = (-t268 * t452 - (pkin(2) * t64 - t22 * t268) * t464) * t118;
t15 = t120 * t354 * t450 + (-t248 * t328 + (-t156 * t368 + (pkin(2) * t272 + t456) * t248) * t66) * t135 * t66;
t14 = t119 * t354 * t451 + (-t248 * t329 + (-t155 * t370 + (pkin(2) * t270 + t457) * t248) * t65) * t134 * t65;
t13 = t118 * t354 * t452 + (-t248 * t330 + (-t154 * t372 + (pkin(2) * t268 + t458) * t248) * t64) * t133 * t64;
t12 = (((t248 * t66 + t60 * t361) * t456 + ((-t334 + t415) * t267 + pkin(2) * t418) * t362 + t24 * t248) * t60 + (t66 * t361 + (t244 * t248 - t346 * t368 - t248) * t60) * t462) * t120;
t11 = (((t248 * t65 + t59 * t363) * t457 + ((-t335 + t416) * t265 + pkin(2) * t419) * t364 + t23 * t248) * t59 + (t65 * t363 + (t243 * t248 - t348 * t370 - t248) * t59) * t463) * t119;
t10 = (((t248 * t64 + t58 * t365) * t458 + ((-t336 + t417) * t263 + pkin(2) * t420) * t366 + t22 * t248) * t58 + (t64 * t365 + (t242 * t248 - t350 * t372 - t248) * t58) * t464) * t118;
t9 = -t144 * t18 - t171 * t12 - Ifges(3,3) * t15 + t57 * (pkin(2) * t445 + t468) + (t474 * mrSges(3,1) + t285 * mrSges(3,2)) * t272 + (t285 * mrSges(3,1) - t474 * mrSges(3,2)) * t266;
t8 = -t143 * t17 - t170 * t11 - Ifges(3,3) * t14 + t56 * (pkin(2) * t446 + t467) + (t473 * mrSges(3,1) + t286 * mrSges(3,2)) * t270 + (t286 * mrSges(3,1) - t473 * mrSges(3,2)) * t264;
t7 = -t142 * t16 - t169 * t10 - Ifges(3,3) * t13 + t55 * (pkin(2) * t447 + t466) + (t472 * mrSges(3,1) + t287 * mrSges(3,2)) * t268 + (t287 * mrSges(3,1) - t472 * mrSges(3,2)) * t262;
t6 = -t18 * t391 - t141 * t12 - t171 * t15 + ((t461 * t266 + t460 * t272) * t66 + (t275 * t266 + t468) * t60) * t478 + (-t183 * t471 - t295 * t212) * t273 + (t295 * t183 - t212 * t471) * t267;
t5 = -t17 * t392 - t140 * t11 - t170 * t14 + ((t461 * t264 + t460 * t270) * t65 + (t275 * t264 + t467) * t59) * t479 + (-t182 * t470 - t297 * t212) * t271 + (t297 * t182 - t212 * t470) * t265;
t4 = -t16 * t393 - t139 * t10 - t169 * t13 + ((t461 * t262 + t460 * t268) * t64 + (t275 * t262 + t466) * t58) * t480 + (-t181 * t469 - t299 * t212) * t269 + (t299 * t181 - t212 * t469) * t263;
t3 = -t144 * t15 - t63 * t248 * t304 + (-t150 * t12 + (-t236 * t57 - t310 * (t57 + t63)) * t267 + (t212 * t60 + t304 * t478) * t418) * t246 + (-t18 - t147) * t240;
t2 = -t143 * t14 - t62 * t248 * t305 + (-t149 * t11 + (-t236 * t56 - t311 * (t56 + t62)) * t265 + (t212 * t59 + t305 * t479) * t419) * t246 + (-t17 - t146) * t240;
t1 = -t142 * t13 - t61 * t248 * t306 + (-t148 * t10 + (-t236 * t55 - t312 * (t55 + t61)) * t263 + (t212 * t58 + t306 * t480) * t420) * t246 + (-t16 - t145) * t240;
t19 = [-g(1) * m(4) + (t75 * t3 - t6 * t421) * t120 + (t74 * t2 - t5 * t422) * t119 + (t73 * t1 - t4 * t423) * t118 + (t9 * t325 + t8 * t326 + t7 * t327) * t283 + (t49 * t441 + t50 * t439 + t51 * t437 + t46 * t434 + t47 * t432 + t48 * t430 + (t54 * t425 + t53 * t427 + t52 * t429) * t283) * t260 + (t49 * t440 + t50 * t438 + t51 * t436 + t46 * t435 + t47 * t433 + t48 * t431 + (t54 * t424 + t53 * t426 + t52 * t428) * t283) * t259 + (m(4) + (-t48 * t421 + t51 * t75) * t120 + (-t47 * t422 + t50 * t74) * t119 + (-t46 * t423 + t49 * t73) * t118 + (t54 * t325 + t53 * t326 + t52 * t327) * t283) * t261; t1 * t441 + t2 * t439 + t3 * t437 + t4 * t434 + t5 * t432 + t6 * t430 - g(2) * m(4) + (t9 * t425 + t8 * t427 + t7 * t429) * t283 + ((-t32 * t421 + t38 * t75) * t120 + (-t30 * t422 + t36 * t74) * t119 + (-t28 * t423 + t34 * t73) * t118 + (t44 * t325 + t42 * t326 + t40 * t327) * t283) * t261 + (t34 * t440 + t36 * t438 + t38 * t436 + t28 * t435 + t30 * t433 + t32 * t431 + (t40 * t428 + t42 * t426 + t44 * t424) * t283) * t259 + (t34 * t441 + t36 * t439 + t38 * t437 + t28 * t434 + t30 * t432 + t32 * t430 + m(4) + (t40 * t429 + t42 * t427 + t44 * t425) * t283) * t260; t1 * t440 + t2 * t438 + t3 * t436 + t4 * t435 + t5 * t433 + t6 * t431 - g(3) * m(4) + (t9 * t424 + t8 * t426 + t7 * t428) * t283 + ((-t33 * t421 + t39 * t75) * t120 + (-t31 * t422 + t37 * t74) * t119 + (-t29 * t423 + t35 * t73) * t118 + (t45 * t325 + t43 * t326 + t41 * t327) * t283) * t261 + (t35 * t441 + t37 * t439 + t39 * t437 + t29 * t434 + t31 * t432 + t33 * t430 + (t41 * t429 + t45 * t425 + t43 * t427) * t283) * t260 + (t35 * t440 + t37 * t438 + t39 * t436 + t29 * t435 + t31 * t433 + t33 * t431 + m(4) + (t41 * t428 + t45 * t424 + t43 * t426) * t283) * t259;];
tauX  = t19;
