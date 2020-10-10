% Calculate vector of inverse dynamics forces for parallel robot
% P3RRRRR7V1G2A0
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d3,d4]';
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
% Datum: 2020-08-07 03:52
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRRRR7V1G2A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V1G2A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRRRR7V1G2A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRRRR7V1G2A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V1G2A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR7V1G2A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRRRR7V1G2A0_invdyn_para_pf_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR7V1G2A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR7V1G2A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRRRR7V1G2A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V1G2A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V1G2A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:42:11
% EndTime: 2020-08-07 03:42:36
% DurationCPUTime: 26.37s
% Computational Cost: add. (147054->1027), mult. (252279->1486), div. (26799->17), fcn. (190029->159), ass. (0->649)
t731 = (-pkin(5) - pkin(4));
t732 = pkin(1) * pkin(2);
t557 = (t731 * t732);
t772 = 4 * t557;
t771 = 2 * t731;
t770 = (m(3) * pkin(4));
t375 = 4 * qJ(3,3);
t380 = 3 * qJ(2,3);
t236 = cos((t375 + t380));
t377 = 2 * qJ(3,3);
t307 = t377 + qJ(2,3);
t239 = cos(t307);
t769 = t236 + t239;
t378 = -2 * qJ(3,3);
t768 = cos((qJ(2,3) + t378)) + cos(t380);
t386 = -2 * qJ(3,2);
t388 = 3 * qJ(2,2);
t767 = cos((t386 + qJ(2,2))) + cos(t388);
t385 = 2 * qJ(3,2);
t315 = t385 + qJ(2,2);
t247 = cos(t315);
t383 = 4 * qJ(3,2);
t250 = cos((t388 + t383));
t766 = t250 + t247;
t393 = 2 * qJ(3,1);
t325 = t393 + qJ(2,1);
t256 = cos(t325);
t391 = 4 * qJ(3,1);
t396 = 3 * qJ(2,1);
t258 = cos((t396 + t391));
t765 = t258 + t256;
t394 = -2 * qJ(3,1);
t764 = cos((qJ(2,1) + t394)) + cos(t396);
t238 = cos((t377 + t380));
t353 = cos(qJ(2,3));
t763 = t353 + t238;
t246 = cos((t385 + t388));
t356 = cos(qJ(2,2));
t762 = t356 + t246;
t259 = cos((t396 + t393));
t359 = cos(qJ(2,1));
t761 = t359 + t259;
t747 = -2 * pkin(2);
t402 = pkin(2) ^ 2;
t406 = 1 / pkin(1);
t590 = t402 * t406;
t760 = 2 * t590;
t759 = 4 * t590;
t364 = xDP(3);
t312 = qJ(2,3) + qJ(3,3);
t221 = sin(t312);
t381 = 2 * qJ(2,3);
t306 = t377 + t381;
t310 = t381 + qJ(3,3);
t313 = qJ(2,3) - qJ(3,3);
t344 = sin(qJ(2,3));
t354 = cos(qJ(1,3));
t382 = -2 * qJ(2,3);
t577 = qJ(3,3) + qJ(1,3);
t578 = -qJ(3,3) + qJ(1,3);
t683 = ((-sin(qJ(2,3) - t578) + sin(qJ(2,3) + t577)) * t771 + (-cos(t382 + t378 + qJ(1,3)) - cos(qJ(1,3) + t306) - 0.2e1 * t354) * pkin(2) + (-cos(t382 + t578) - cos(qJ(1,3) + t310) - cos(t578) - cos(t577)) * pkin(1)) / ((-sin(t307) + t344) * pkin(2) + (sin(t313) - t221) * pkin(1));
t484 = t683 / 0.2e1;
t437 = t406 * t484;
t365 = xDP(2);
t352 = cos(qJ(3,3));
t721 = pkin(2) * t352;
t200 = pkin(1) + t721;
t343 = sin(qJ(3,3));
t611 = t343 * t344;
t421 = pkin(2) * t611 - t200 * t353;
t145 = 0.1e1 / t421;
t296 = 0.1e1 / t343;
t622 = t296 * t406;
t509 = t145 * t622;
t299 = t352 ^ 2;
t345 = sin(qJ(1,3));
t596 = t354 * t731;
t142 = t596 * t611 + (t299 - 0.1e1) * t345 * pkin(2);
t337 = legFrame(3,2);
t268 = sin(t337);
t271 = cos(t337);
t300 = t353 ^ 2;
t494 = t345 * t611;
t746 = 2 * pkin(2);
t418 = pkin(1) * t494 + (t494 * t746 - t596) * t352;
t610 = t343 * t352;
t542 = pkin(2) * t610;
t562 = 0.2e1 * t721;
t632 = (pkin(1) + t562) * t343;
t276 = pkin(1) * t352;
t175 = t299 * t746 - pkin(2) + t276;
t637 = t175 * t345;
t638 = t175 * t344;
t82 = (-t268 * t637 + t271 * t632) * t300 + (t268 * t418 + t271 * t638) * t353 + t142 * t268 - t271 * t542;
t452 = t82 * t509;
t76 = t365 * t452;
t366 = xDP(1);
t83 = (t268 * t632 + t271 * t637) * t300 + (t268 * t638 - t271 * t418) * t353 - t142 * t271 - t268 * t542;
t451 = t83 * t509;
t77 = t366 * t451;
t714 = -t76 - t77;
t49 = t364 * t437 + t714;
t730 = pkin(1) * t49;
t321 = qJ(2,2) + qJ(3,2);
t228 = sin(t321);
t389 = 2 * qJ(2,2);
t317 = t389 + t385;
t318 = t389 + qJ(3,2);
t322 = qJ(2,2) - qJ(3,2);
t347 = sin(qJ(2,2));
t357 = cos(qJ(1,2));
t390 = -2 * qJ(2,2);
t579 = qJ(3,2) + qJ(1,2);
t580 = -qJ(3,2) + qJ(1,2);
t682 = ((-sin(qJ(2,2) - t580) + sin(qJ(2,2) + t579)) * t771 + (-cos(t390 + t386 + qJ(1,2)) - cos(qJ(1,2) + t317) - 0.2e1 * t357) * pkin(2) + (-cos(t390 + t580) - cos(qJ(1,2) + t318) - cos(t580) - cos(t579)) * pkin(1)) / ((-sin(t315) + t347) * pkin(2) + (sin(t322) - t228) * pkin(1));
t483 = t682 / 0.2e1;
t436 = t406 * t483;
t355 = cos(qJ(3,2));
t720 = pkin(2) * t355;
t202 = pkin(1) + t720;
t346 = sin(qJ(3,2));
t606 = t346 * t347;
t420 = pkin(2) * t606 - t202 * t356;
t146 = 0.1e1 / t420;
t297 = 0.1e1 / t346;
t620 = t297 * t406;
t508 = t146 * t620;
t301 = t355 ^ 2;
t348 = sin(qJ(1,2));
t595 = t357 * t731;
t143 = t595 * t606 + (t301 - 0.1e1) * t348 * pkin(2);
t338 = legFrame(2,2);
t269 = sin(t338);
t272 = cos(t338);
t302 = t356 ^ 2;
t493 = t348 * t606;
t417 = pkin(1) * t493 + (t493 * t746 - t595) * t355;
t605 = t346 * t355;
t540 = pkin(2) * t605;
t561 = 0.2e1 * t720;
t630 = (pkin(1) + t561) * t346;
t277 = pkin(1) * t355;
t176 = t301 * t746 - pkin(2) + t277;
t635 = t176 * t348;
t636 = t176 * t347;
t84 = (-t269 * t635 + t272 * t630) * t302 + (t269 * t417 + t272 * t636) * t356 + t143 * t269 - t272 * t540;
t450 = t84 * t508;
t78 = t365 * t450;
t85 = (t269 * t630 + t272 * t635) * t302 + (t269 * t636 - t272 * t417) * t356 - t143 * t272 - t269 * t540;
t449 = t85 * t508;
t79 = t366 * t449;
t713 = -t78 - t79;
t50 = t364 * t436 + t713;
t729 = pkin(1) * t50;
t330 = qJ(2,1) + qJ(3,1);
t235 = sin(t330);
t397 = 2 * qJ(2,1);
t324 = t393 + t397;
t328 = t397 + qJ(3,1);
t331 = qJ(2,1) - qJ(3,1);
t350 = sin(qJ(2,1));
t360 = cos(qJ(1,1));
t398 = -2 * qJ(2,1);
t581 = qJ(3,1) + qJ(1,1);
t582 = -qJ(3,1) + qJ(1,1);
t681 = ((-sin(qJ(2,1) - t582) + sin(qJ(2,1) + t581)) * t771 + (-cos(t398 + t394 + qJ(1,1)) - cos(qJ(1,1) + t324) - 0.2e1 * t360) * pkin(2) + (-cos(t398 + t582) - cos(qJ(1,1) + t328) - cos(t582) - cos(t581)) * pkin(1)) / ((-sin(t325) + t350) * pkin(2) + (sin(t331) - t235) * pkin(1));
t482 = t681 / 0.2e1;
t435 = t406 * t482;
t358 = cos(qJ(3,1));
t719 = pkin(2) * t358;
t204 = pkin(1) + t719;
t349 = sin(qJ(3,1));
t601 = t349 * t350;
t419 = pkin(2) * t601 - t204 * t359;
t147 = 0.1e1 / t419;
t298 = 0.1e1 / t349;
t618 = t298 * t406;
t507 = t147 * t618;
t303 = t358 ^ 2;
t351 = sin(qJ(1,1));
t594 = t360 * t731;
t144 = t594 * t601 + (t303 - 0.1e1) * t351 * pkin(2);
t339 = legFrame(1,2);
t270 = sin(t339);
t273 = cos(t339);
t304 = t359 ^ 2;
t492 = t351 * t601;
t416 = pkin(1) * t492 + (t492 * t746 - t594) * t358;
t600 = t349 * t358;
t538 = pkin(2) * t600;
t560 = 0.2e1 * t719;
t628 = (pkin(1) + t560) * t349;
t278 = pkin(1) * t358;
t177 = t303 * t746 - pkin(2) + t278;
t633 = t177 * t351;
t634 = t177 * t350;
t86 = (-t270 * t633 + t273 * t628) * t304 + (t270 * t416 + t273 * t634) * t359 + t144 * t270 - t273 * t538;
t448 = t86 * t507;
t80 = t365 * t448;
t87 = (t270 * t628 + t273 * t633) * t304 + (t270 * t634 - t273 * t416) * t359 - t144 * t273 - t270 * t538;
t447 = t87 * t507;
t81 = t366 * t447;
t712 = -t80 - t81;
t51 = t364 * t435 + t712;
t728 = pkin(1) * t51;
t403 = 1 / pkin(2);
t623 = t296 * t403;
t136 = t345 * t731 + t354 * t421;
t656 = t136 * t364;
t133 = -t345 * t421 + t596;
t609 = t343 * t353;
t541 = pkin(2) * t609;
t631 = t200 * t344;
t148 = t541 + t631;
t116 = -t133 * t271 - t148 * t268;
t672 = t116 * t366;
t115 = t133 * t268 - t148 * t271;
t674 = t115 * t365;
t411 = (t656 + t672 + t674) * t623;
t73 = t406 * t411;
t34 = t73 + t49;
t621 = t297 * t403;
t137 = t348 * t731 + t357 * t420;
t655 = t137 * t364;
t134 = -t348 * t420 + t595;
t604 = t346 * t356;
t539 = pkin(2) * t604;
t629 = t202 * t347;
t149 = t539 + t629;
t118 = -t134 * t272 - t149 * t269;
t668 = t118 * t366;
t117 = t134 * t269 - t149 * t272;
t670 = t117 * t365;
t410 = (t655 + t668 + t670) * t621;
t74 = t406 * t410;
t35 = t74 + t50;
t619 = t298 * t403;
t138 = t351 * t731 + t360 * t419;
t654 = t138 * t364;
t135 = -t351 * t419 + t594;
t599 = t349 * t359;
t537 = pkin(2) * t599;
t627 = t204 * t350;
t150 = t537 + t627;
t120 = -t135 * t273 - t150 * t270;
t664 = t120 * t366;
t119 = t135 * t270 - t150 * t273;
t666 = t119 * t365;
t409 = (t654 + t664 + t666) * t619;
t75 = t406 * t409;
t36 = t75 + t51;
t589 = t403 * t406;
t498 = t298 * t589;
t757 = (t664 / 0.2e1 + t666 / 0.2e1 + t654 / 0.2e1) * t498;
t499 = t297 * t589;
t756 = (t668 / 0.2e1 + t670 / 0.2e1 + t655 / 0.2e1) * t499;
t500 = t296 * t589;
t755 = (t672 / 0.2e1 + t674 / 0.2e1 + t656 / 0.2e1) * t500;
t266 = pkin(4) * mrSges(3,2) - Ifges(3,6);
t267 = mrSges(3,1) * pkin(4) - Ifges(3,5);
t576 = t343 * t266 - t267 * t352;
t575 = t346 * t266 - t267 * t355;
t574 = t349 * t266 - t267 * t358;
t711 = mrSges(3,1) * t343;
t207 = pkin(1) * t711;
t335 = Ifges(3,4) - Ifges(2,4);
t754 = t335 + t207;
t710 = mrSges(3,1) * t346;
t209 = pkin(1) * t710;
t753 = t335 + t209;
t709 = mrSges(3,1) * t349;
t211 = pkin(1) * t709;
t752 = t335 + t211;
t522 = mrSges(3,3) + t770;
t193 = pkin(1) * t522 - Ifges(2,5);
t750 = -0.3e1 * pkin(1);
t749 = -0.2e1 * pkin(1);
t748 = -6 * pkin(2);
t745 = -2 * Ifges(3,4);
t25 = t49 + t755;
t744 = 0.4e1 * t25;
t26 = t50 + t756;
t743 = 0.4e1 * t26;
t27 = t51 + t757;
t742 = 0.4e1 * t27;
t242 = cos(t312);
t169 = 0.1e1 / (pkin(1) * t353 + pkin(2) * t242);
t112 = (-t345 * t364 + (-t268 * t365 + t271 * t366) * t354) * t169;
t741 = 0.2e1 * t112;
t252 = cos(t321);
t170 = 0.1e1 / (pkin(1) * t356 + pkin(2) * t252);
t113 = (-t348 * t364 + (-t269 * t365 + t272 * t366) * t357) * t170;
t740 = 0.2e1 * t113;
t260 = cos(t330);
t171 = 0.1e1 / (pkin(1) * t359 + pkin(2) * t260);
t114 = (-t351 * t364 + (-t270 * t365 + t273 * t366) * t360) * t171;
t739 = 0.2e1 * t114;
t334 = Ifges(3,2) - Ifges(3,1);
t405 = pkin(1) ^ 2;
t373 = m(3) * t405;
t194 = -Ifges(2,1) + Ifges(2,2) + t373 - t334;
t738 = -0.2e1 * t194;
t705 = mrSges(3,2) * t343;
t206 = pkin(1) * t705;
t737 = -0.2e1 * t206;
t704 = mrSges(3,2) * t346;
t208 = pkin(1) * t704;
t736 = -0.2e1 * t208;
t703 = mrSges(3,2) * t349;
t210 = pkin(1) * t703;
t735 = -0.2e1 * t210;
t733 = -pkin(1) / 0.2e1;
t362 = pkin(1) * mrSges(3,1);
t726 = -t402 / 0.2e1;
t725 = t402 / 0.2e1;
t724 = -t405 / 0.2e1;
t723 = pkin(1) * (-(3 * t402) + t405);
t722 = pkin(1) * t402;
t718 = pkin(2) * t405;
t28 = t402 * t34;
t43 = t405 * t49;
t717 = t28 + t43;
t29 = t402 * t35;
t44 = t405 * t50;
t716 = t29 + t44;
t30 = t402 * t36;
t45 = t405 * t51;
t715 = t30 + t45;
t708 = mrSges(3,1) * t352;
t707 = mrSges(3,1) * t355;
t706 = mrSges(3,1) * t358;
t702 = Ifges(3,4) * t299;
t701 = Ifges(3,4) * t301;
t700 = Ifges(3,4) * t303;
t699 = Ifges(3,4) * t343;
t698 = Ifges(3,4) * t346;
t697 = Ifges(3,4) * t349;
t696 = t145 * t82;
t695 = t145 * t83;
t694 = t146 * t84;
t693 = t146 * t85;
t692 = t147 * t86;
t691 = t147 * t87;
t690 = t34 * t299;
t689 = t35 * t301;
t688 = t36 * t303;
t684 = m(3) * pkin(1) + mrSges(2,1);
t680 = t112 * t300;
t679 = t112 * t731;
t678 = t113 * t302;
t677 = t113 * t731;
t676 = t114 * t304;
t675 = t114 * t731;
t673 = t115 * t403;
t671 = t116 * t403;
t669 = t117 * t403;
t667 = t118 * t403;
t665 = t119 * t403;
t663 = t120 * t403;
t626 = t267 * t343;
t458 = -t266 * t352 - t626;
t121 = (-t193 + t576) * t344 + (Ifges(2,6) + t458) * t353;
t662 = t121 * t145;
t625 = t267 * t346;
t457 = -t266 * t355 - t625;
t122 = (-t193 + t575) * t347 + (Ifges(2,6) + t457) * t356;
t661 = t122 * t146;
t624 = t267 * t349;
t456 = -t266 * t358 - t624;
t123 = (-t193 + t574) * t350 + (Ifges(2,6) + t456) * t359;
t660 = t123 * t147;
t124 = t344 * t576 + t353 * t458;
t659 = t124 * t403;
t125 = t347 * t575 + t356 * t457;
t658 = t125 * t403;
t126 = t350 * t574 + t359 * t456;
t657 = t126 * t403;
t526 = Ifges(2,3) + Ifges(3,3) + t373;
t550 = pkin(1) * t708;
t154 = t737 + t526 + 0.2e1 * t550;
t653 = t145 * t154;
t160 = Ifges(3,3) - t206 + t550;
t652 = t145 * t160;
t549 = pkin(1) * t707;
t155 = t736 + t526 + 0.2e1 * t549;
t651 = t146 * t155;
t161 = Ifges(3,3) - t208 + t549;
t650 = t146 * t161;
t548 = pkin(1) * t706;
t156 = t735 + t526 + 0.2e1 * t548;
t649 = t147 * t156;
t162 = Ifges(3,3) - t210 + t548;
t648 = t147 * t162;
t647 = t160 * t403;
t646 = t161 * t403;
t645 = t162 * t403;
t644 = t169 * t345;
t643 = t169 * t354;
t642 = t170 * t348;
t641 = t170 * t357;
t640 = t171 * t351;
t639 = t171 * t360;
t617 = t299 * t402;
t616 = t301 * t402;
t615 = t303 * t402;
t614 = t334 * t343;
t613 = t334 * t346;
t612 = t334 * t349;
t608 = t344 * t353;
t607 = t344 * t731;
t603 = t347 * t356;
t602 = t347 * t731;
t598 = t350 * t359;
t597 = t350 * t731;
t593 = t731 * t402;
t592 = t731 * t405;
t591 = (-pkin(2) + pkin(1)) * (pkin(1) + pkin(2));
t308 = t378 + t381;
t404 = pkin(1) * t405;
t588 = t404 * sin(t308);
t320 = t389 + t386;
t587 = t404 * sin(t320);
t326 = t394 + t397;
t586 = t404 * sin(t326);
t311 = t381 - qJ(3,3);
t585 = t405 * sin(t311);
t319 = t389 - qJ(3,2);
t584 = t405 * sin(t319);
t329 = t397 - qJ(3,1);
t583 = t405 * sin(t329);
t332 = t731 ^ 2;
t573 = 0.4e1 / 0.3e1 * t332 + t405;
t572 = 0.2e1 / 0.3e1 * t402 + t405;
t275 = -Ifges(3,4) / 0.2e1 + Ifges(2,4) / 0.2e1;
t280 = t402 + t405;
t571 = t405 - t402;
t570 = 0.2e1 * pkin(1) * (0.9e1 / 0.2e1 * t402 + t405 + (2 * t332));
t569 = ((3 * t402) + t573) * t750;
t568 = ((2 * t402) + t573) * t750;
t567 = -0.3e1 * t722;
t566 = 0.5e1 * t722;
t565 = -0.9e1 * pkin(2) * (0.8e1 / 0.9e1 * t332 + t572);
t564 = (0.2e1 / 0.3e1 * t332 + t572) * t748;
t205 = t332 + t402;
t563 = 4 * pkin(2) * t205;
t559 = 0.3e1 * t718;
t558 = t406 * t747;
t553 = -0.4e1 * t690;
t552 = -0.4e1 * t689;
t551 = -0.4e1 * t688;
t547 = 2 * t593;
t546 = -0.2e1 * t592;
t536 = pkin(2) * t585;
t535 = pkin(2) * t584;
t534 = pkin(2) * t583;
t376 = 3 * qJ(3,3);
t533 = sin(t376) * t718;
t384 = 3 * qJ(3,2);
t532 = sin(t384) * t718;
t392 = 3 * qJ(3,1);
t531 = sin(t392) * t718;
t530 = t362 / 0.2e1;
t529 = t25 * t733;
t528 = t26 * t733;
t527 = t27 * t733;
t525 = t34 * t614;
t524 = t35 * t613;
t523 = t36 * t612;
t521 = t364 * t683;
t520 = t364 * t682;
t519 = t364 * t681;
t512 = t136 * t623;
t511 = t137 * t621;
t510 = t138 * t619;
t506 = t268 * t643;
t505 = t271 * t643;
t504 = t269 * t641;
t503 = t272 * t641;
t502 = t270 * t639;
t501 = t273 * t639;
t497 = t334 * t610;
t496 = t334 * t605;
t495 = t334 * t600;
t491 = t731 * t591;
t489 = 2 * t557;
t488 = -4 * t557;
t479 = 0.2e1 / 0.3e1 * t589;
t432 = t366 * t479;
t433 = t365 * t479;
t434 = t364 * t479;
t487 = (t115 * t433 + t116 * t432 + t136 * t434) * t296;
t486 = (t117 * t433 + t118 * t432 + t137 * t434) * t297;
t485 = (t119 * t433 + t120 * t432 + t138 * t434) * t298;
t481 = t364 * t406 / 0.3e1;
t480 = t591 / 0.2e1;
t478 = t402 * t714;
t477 = t402 * t713;
t476 = t402 * t712;
t469 = t591 * t747;
t401 = pkin(2) * t402;
t468 = t401 - 0.3e1 * t718;
t467 = t112 * t557;
t466 = t113 * t557;
t465 = t114 * t557;
t305 = t376 + t381;
t215 = sin(t305);
t216 = sin(t306);
t309 = t381 + t375;
t218 = sin(t309);
t219 = sin(t310);
t237 = cos((t376 + qJ(2,3)));
t240 = cos((t380 + t376));
t241 = cos((t380 + qJ(3,3)));
t243 = cos(t313);
t263 = t724 + t402;
t264 = t405 + t726;
t283 = sin(t377);
t284 = sin(t381);
t379 = 4 * qJ(2,3);
t10 = (((-t585 / 0.2e1 + t215 * t725 + t219 * t480 + t280 * t343 + (t216 + t283 - t284) * t732) * t34 - ((-t238 / 0.4e1 + t236 / 0.4e1 + t239 / 0.4e1 - t353 / 0.4e1) * pkin(2) + (t240 / 0.4e1 - t241 / 0.4e1 + t237 / 0.4e1 - t242 / 0.4e1) * pkin(1)) * t679) * t73 * t558 + ((t73 * t237 * t489 + t763 * (t43 - t28) * t771 + t769 * t34 * t547 + ((t218 + sin((t375 + t379))) * t567 + t468 * sin((t379 + t376)) + t404 * sin(t379) - t401 * sin((5 * qJ(3,3) + t379)) + sin((t379 + qJ(3,3))) * t559 - sin((t377 + t379)) * t723 + t215 * t564 + t216 * t568 + t219 * t563 + t283 * t569 + t284 * t570 + t343 * t565 - 0.3e1 * t533 + 0.6e1 * t536 + t588) * t112 + (t243 * t488 + t768 * t546) * t49 + ((t49 - t755) * t242 - (-t240 + t241) * t25) * t772) * t112 / 0.4e1 + (-pkin(1) * (t411 * t759 + t521 * t760 + t43 + 0.4e1 * t478) * t216 + t717 * t215 * t747 + t536 * t744 - t49 * t218 * t722 + t49 * t588 + t34 * t219 * t469 + (0.3e1 * t478 + t43 + (0.3e1 / 0.2e1 * t521 + 0.2e1 * t411) * t590) * t283 * t749 + ((0.4e1 / 0.5e1 * t672 + 0.4e1 / 0.5e1 * t674 + 0.4e1 / 0.5e1 * t656) * t500 + t49) * t284 * t566 - 0.2e1 * t49 * t533 + t343 * ((-0.2e1 / 0.3e1 * t77 - 0.2e1 / 0.3e1 * t76 + t481 * t683 + t487) * t402 + t405 * (t487 + t49)) * t748 + 0.2e1 * (t240 + t242) * t467 - 0.2e1 * (t243 + t241) * t467 + (t763 * t491 - t768 * t592 + t769 * t593) * t112) * t49 / 0.2e1) * t406) / (t263 * cos(t306) + cos(t309) * t726 + cos(t308) * t724 + t264 * cos(t381) - t280 * cos(t377) + (-cos(t311) - cos(t305) + 0.2e1 * cos(t310) - cos(t376) + t352) * t732 + t280);
t109 = t112 ^ 2;
t13 = (t221 * t679 + t34 * (t276 + pkin(2))) * t73 * t622 + (((0.2e1 * t25 * t276 - (-t344 * t352 - t609) * t679) * pkin(2) + t717) * t49 + ((pkin(1) * t562 + t571 + 0.2e1 * t617) * t680 + 0.2e1 * t607 * t730 + (-0.2e1 * t541 * t631 + t205 - t617) * t112 + (t541 * t731 + t607 * t721) * t34) * t112) * t500;
t16 = (t34 * t221 * t747 - 0.2e1 * t344 * t730 - t679) * t112 * t169;
t163 = g(1) * t268 + g(2) * t271;
t336 = 0.2e1 * t334;
t196 = t336 * t299;
t212 = 0.4e1 * t699;
t361 = pkin(1) * mrSges(3,2);
t399 = 2 * Ifges(3,4);
t166 = g(1) * t271 - g(2) * t268;
t428 = g(3) * t354 + t166 * t345;
t46 = t49 ^ 2;
t97 = 0.2e1 * t109 * t702;
t455 = (-Ifges(3,3) * t13 - t160 * t10 - t124 * t16 + t97 + (t221 * t163 + t242 * t428 + t276 * t46) * mrSges(3,2) + (t343 * pkin(1) * t46 - t163 * t242 + t221 * t428) * mrSges(3,1) + ((-0.4e1 * t702 + (t336 * t343 + t361) * t352 + t207 + t399) * t300 + (t196 + (t212 + t362) * t352 - t206 - t334) * t608 - t497 - Ifges(3,4)) * t109) * t500;
t110 = t113 ^ 2;
t314 = t384 + t389;
t222 = sin(t314);
t316 = t389 + t383;
t223 = sin(t316);
t224 = sin(t317);
t225 = sin(t318);
t245 = cos((t384 + qJ(2,2)));
t248 = cos((qJ(3,2) + t388));
t251 = cos((t388 + t384));
t253 = cos(t322);
t286 = sin(t385);
t287 = sin(t389);
t387 = 4 * qJ(2,2);
t12 = (((-t584 / 0.2e1 + t222 * t725 + t225 * t480 + t280 * t346 + (t224 + t286 - t287) * t732) * t35 - ((t250 / 0.4e1 - t246 / 0.4e1 + t247 / 0.4e1 - t356 / 0.4e1) * pkin(2) + (t251 / 0.4e1 - t248 / 0.4e1 + t245 / 0.4e1 - t252 / 0.4e1) * pkin(1)) * t677) * t74 * t558 + ((t74 * t245 * t489 + t762 * (t44 - t29) * t771 + t766 * t35 * t547 + ((t223 + sin((t387 + t383))) * t567 + t468 * sin((t384 + t387)) + t404 * sin(t387) - t401 * sin((5 * qJ(3,2) + t387)) + sin((t387 + qJ(3,2))) * t559 - sin((t387 + t385)) * t723 + t222 * t564 + t224 * t568 + t225 * t563 + t286 * t569 + t287 * t570 + t346 * t565 - 0.3e1 * t532 + 0.6e1 * t535 + t587) * t113 + (t253 * t488 + t767 * t546) * t50 + ((t50 - t756) * t252 - (-t251 + t248) * t26) * t772) * t113 / 0.4e1 + (-pkin(1) * (t410 * t759 + t520 * t760 + t44 + 0.4e1 * t477) * t224 + t716 * t222 * t747 + t535 * t743 + t50 * t587 - t50 * t223 * t722 + t35 * t225 * t469 + (0.3e1 * t477 + t44 + (0.3e1 / 0.2e1 * t520 + 0.2e1 * t410) * t590) * t286 * t749 + ((0.4e1 / 0.5e1 * t668 + 0.4e1 / 0.5e1 * t670 + 0.4e1 / 0.5e1 * t655) * t499 + t50) * t287 * t566 - 0.2e1 * t50 * t532 + t346 * ((-0.2e1 / 0.3e1 * t79 - 0.2e1 / 0.3e1 * t78 + t481 * t682 + t486) * t402 + t405 * (t486 + t50)) * t748 + 0.2e1 * (t251 + t252) * t466 - 0.2e1 * (t253 + t248) * t466 + (t762 * t491 - t767 * t592 + t766 * t593) * t113) * t50 / 0.2e1) * t406) / (t263 * cos(t317) + cos(t320) * t724 + cos(t316) * t726 + t264 * cos(t389) - t280 * cos(t385) + (-cos(t319) - cos(t314) + 0.2e1 * cos(t318) - cos(t384) + t355) * t732 + t280);
t14 = (t228 * t677 + t35 * (t277 + pkin(2))) * t74 * t620 + (((0.2e1 * t26 * t277 - (-t347 * t355 - t604) * t677) * pkin(2) + t716) * t50 + ((pkin(1) * t561 + t571 + 0.2e1 * t616) * t678 + 0.2e1 * t602 * t729 + (-0.2e1 * t539 * t629 + t205 - t616) * t113 + (t539 * t731 + t602 * t720) * t35) * t113) * t499;
t164 = g(1) * t269 + g(2) * t272;
t17 = (t35 * t228 * t747 - 0.2e1 * t347 * t729 - t677) * t113 * t170;
t197 = t336 * t301;
t213 = 0.4e1 * t698;
t167 = g(1) * t272 - g(2) * t269;
t427 = g(3) * t357 + t167 * t348;
t47 = t50 ^ 2;
t98 = 0.2e1 * t110 * t701;
t454 = (-Ifges(3,3) * t14 - t161 * t12 - t125 * t17 + t98 + (t228 * t164 + t252 * t427 + t277 * t47) * mrSges(3,2) + (t346 * pkin(1) * t47 - t164 * t252 + t228 * t427) * mrSges(3,1) + ((-0.4e1 * t701 + (t336 * t346 + t361) * t355 + t209 + t399) * t302 + (t197 + (t213 + t362) * t355 - t208 - t334) * t603 - t496 - Ifges(3,4)) * t110) * t499;
t323 = t392 + t397;
t229 = sin(t323);
t230 = sin(t324);
t327 = t397 + t391;
t232 = sin(t327);
t233 = sin(t328);
t254 = cos((t392 + t396));
t255 = cos((t392 + qJ(2,1)));
t257 = cos((qJ(3,1) + t396));
t261 = cos(t331);
t289 = sin(t393);
t290 = sin(t397);
t395 = 4 * qJ(2,1);
t11 = (((-t583 / 0.2e1 + t229 * t725 + t233 * t480 + t280 * t349 + (t230 + t289 - t290) * t732) * t36 - ((-t259 / 0.4e1 + t258 / 0.4e1 + t256 / 0.4e1 - t359 / 0.4e1) * pkin(2) + (t254 / 0.4e1 - t257 / 0.4e1 + t255 / 0.4e1 - t260 / 0.4e1) * pkin(1)) * t675) * t75 * t558 + ((t75 * t255 * t489 + t761 * (t45 - t30) * t771 + t765 * t36 * t547 + ((t232 + sin((t395 + t391))) * t567 + t468 * sin((t392 + t395)) + t404 * sin(t395) - t401 * sin((t395 + 5 * qJ(3,1))) + sin((t395 + qJ(3,1))) * t559 - sin((t395 + t393)) * t723 + t229 * t564 + t230 * t568 + t233 * t563 + t289 * t569 + t290 * t570 + t349 * t565 - 0.3e1 * t531 + 0.6e1 * t534 + t586) * t114 + (t261 * t488 + t764 * t546) * t51 + ((t51 - t757) * t260 - (-t254 + t257) * t27) * t772) * t114 / 0.4e1 + (-pkin(1) * (t409 * t759 + t519 * t760 + t45 + 0.4e1 * t476) * t230 + t715 * t229 * t747 + t534 * t742 - t51 * t232 * t722 + t51 * t586 + t36 * t233 * t469 + (0.3e1 * t476 + t45 + (0.3e1 / 0.2e1 * t519 + 0.2e1 * t409) * t590) * t289 * t749 + ((0.4e1 / 0.5e1 * t664 + 0.4e1 / 0.5e1 * t666 + 0.4e1 / 0.5e1 * t654) * t498 + t51) * t290 * t566 - 0.2e1 * t51 * t531 + t349 * ((-0.2e1 / 0.3e1 * t81 - 0.2e1 / 0.3e1 * t80 + t481 * t681 + t485) * t402 + t405 * (t485 + t51)) * t748 + 0.2e1 * (t254 + t260) * t465 - 0.2e1 * (t261 + t257) * t465 + (t761 * t491 - t764 * t592 + t765 * t593) * t114) * t51 / 0.2e1) * t406) / (t263 * cos(t324) + cos(t327) * t726 + cos(t326) * t724 + t264 * cos(t397) - t280 * cos(t393) + (-cos(t329) - cos(t323) + 0.2e1 * cos(t328) - cos(t392) + t358) * t732 + t280);
t111 = t114 ^ 2;
t15 = (t235 * t675 + t36 * (t278 + pkin(2))) * t75 * t618 + (((0.2e1 * t27 * t278 - (-t350 * t358 - t599) * t675) * pkin(2) + t715) * t51 + ((pkin(1) * t560 + t571 + 0.2e1 * t615) * t676 + 0.2e1 * t597 * t728 + (-0.2e1 * t537 * t627 + t205 - t615) * t114 + (t537 * t731 + t597 * t719) * t36) * t114) * t498;
t165 = g(1) * t270 + g(2) * t273;
t18 = (t36 * t235 * t747 - 0.2e1 * t350 * t728 - t675) * t114 * t171;
t198 = t336 * t303;
t214 = 0.4e1 * t697;
t168 = g(1) * t273 - g(2) * t270;
t426 = g(3) * t360 + t168 * t351;
t48 = t51 ^ 2;
t99 = 0.2e1 * t111 * t700;
t453 = (-Ifges(3,3) * t15 - t162 * t11 - t126 * t18 + t99 + (t235 * t165 + t260 * t426 + t278 * t48) * mrSges(3,2) + (t349 * pkin(1) * t48 - t165 * t260 + t235 * t426) * mrSges(3,1) + ((-0.4e1 * t700 + (t336 * t349 + t361) * t358 + t211 + t399) * t304 + (t198 + (t214 + t362) * t358 - t210 - t334) * t598 - t495 - Ifges(3,4)) * t111) * t498;
t431 = -mrSges(3,2) * t352 - t711;
t430 = -mrSges(3,2) * t355 - t710;
t429 = -mrSges(3,2) * t358 - t709;
t157 = t684 - t705 + t708;
t172 = mrSges(2,2) - t431;
t425 = t157 * t353 - t172 * t344;
t158 = t684 - t704 + t707;
t173 = mrSges(2,2) - t430;
t424 = t158 * t356 - t173 * t347;
t159 = t684 - t703 + t706;
t174 = mrSges(2,2) - t429;
t423 = t159 * t359 - t174 * t350;
t415 = Ifges(2,1) + Ifges(3,2) + Ifges(1,3) + ((2 * mrSges(3,3) + t770) * pkin(4));
t363 = mrSges(1,1) * g(3);
t342 = xDDP(1);
t341 = xDDP(2);
t340 = xDDP(3);
t294 = 0.2e1 * t362;
t279 = -t361 / 0.2e1;
t274 = -Ifges(3,2) / 0.2e1 + Ifges(3,1) / 0.2e1;
t195 = -mrSges(1,2) + mrSges(2,3) + t522;
t180 = t195 * g(3);
t132 = t198 + (t214 + t294) * t358 + t735 + t194;
t131 = t197 + (t213 + t294) * t355 + t736 + t194;
t130 = t196 + (t212 + t294) * t352 + t737 + t194;
t93 = t132 * t304 + 0.4e1 * (t700 + (t274 * t349 + t279) * t358 - t211 / 0.2e1 + t275) * t598 - t334 * t303 + t600 * t745 + t415;
t92 = t131 * t302 + 0.4e1 * (t701 + (t274 * t346 + t279) * t355 - t209 / 0.2e1 + t275) * t603 - t334 * t301 + t605 * t745 + t415;
t91 = t130 * t300 + 0.4e1 * (t702 + (t274 * t343 + t279) * t352 - t207 / 0.2e1 + t275) * t608 - t334 * t299 + t610 * t745 + t415;
t72 = -t126 * t640 + (Ifges(3,3) * t510 + t162 * t482) * t406;
t71 = -t125 * t642 + (Ifges(3,3) * t511 + t161 * t483) * t406;
t70 = -t124 * t644 + (Ifges(3,3) * t512 + t160 * t484) * t406;
t69 = -t123 * t640 + (t156 * t482 + t162 * t510) * t406;
t68 = -t122 * t642 + (t155 * t483 + t161 * t511) * t406;
t67 = -t121 * t644 + (t154 * t484 + t160 * t512) * t406;
t66 = t126 * t501 + (Ifges(3,3) * t663 - t648 * t87) * t618;
t65 = t125 * t503 + (Ifges(3,3) * t667 - t650 * t85) * t620;
t64 = t124 * t505 + (Ifges(3,3) * t671 - t652 * t83) * t622;
t63 = -t126 * t502 + (Ifges(3,3) * t665 - t648 * t86) * t618;
t62 = -t125 * t504 + (Ifges(3,3) * t669 - t650 * t84) * t620;
t61 = -t124 * t506 + (Ifges(3,3) * t673 - t652 * t82) * t622;
t60 = t123 * t501 + (t120 * t645 - t649 * t87) * t618;
t59 = t122 * t503 + (t118 * t646 - t651 * t85) * t620;
t58 = t121 * t505 + (t116 * t647 - t653 * t83) * t622;
t57 = -t123 * t502 + (t119 * t645 - t649 * t86) * t618;
t56 = -t122 * t504 + (t117 * t646 - t651 * t84) * t620;
t55 = -t121 * t506 + (t115 * t647 - t653 * t82) * t622;
t54 = -t93 * t640 + (t123 * t482 + t126 * t510) * t406;
t53 = -t92 * t642 + (t122 * t483 + t125 * t511) * t406;
t52 = -t91 * t644 + (t121 * t484 + t124 * t512) * t406;
t42 = t93 * t501 + (t120 * t657 - t660 * t87) * t618;
t41 = t92 * t503 + (t118 * t658 - t661 * t85) * t620;
t40 = t91 * t505 + (t116 * t659 - t662 * t83) * t622;
t39 = -t93 * t502 + (t119 * t657 - t660 * t86) * t618;
t38 = -t92 * t504 + (t117 * t658 - t661 * t84) * t620;
t37 = -t91 * t506 + (t115 * t659 - t662 * t82) * t622;
t33 = t36 ^ 2;
t32 = t35 ^ 2;
t31 = t34 ^ 2;
t6 = -t123 * t18 - t156 * t11 - t162 * t15 + t99 + (t159 * t426 + t165 * t174) * t350 - t359 * (t159 * t165 - t174 * t426) + t429 * (t111 + (0.2e1 * t51 + t75) * t75) * pkin(1) + (0.2e1 * (-0.2e1 * t700 + (t361 + t612) * t358 + t752) * t304 + t132 * t598 - t495 - t335) * t111;
t5 = -t122 * t17 - t155 * t12 - t161 * t14 + t98 + (t158 * t427 + t164 * t173) * t347 - t356 * (t158 * t164 - t173 * t427) + t430 * (t110 + (0.2e1 * t50 + t74) * t74) * pkin(1) + (0.2e1 * (-0.2e1 * t701 + (t361 + t613) * t355 + t753) * t302 + t131 * t603 - t496 - t335) * t110;
t4 = -t121 * t16 - t154 * t10 - t160 * t13 + t97 + (t157 * t428 + t163 * t172) * t344 - t353 * (t157 * t163 - t172 * t428) + t431 * (t109 + (0.2e1 * t49 + t73) * t73) * pkin(1) + (0.2e1 * (-0.2e1 * t702 + (t361 + t614) * t352 + t754) * t300 + t130 * t608 - t497 - t335) * t109;
t3 = -t93 * t18 - t123 * t11 - t126 * t15 + 0.8e1 * ((-t523 / 0.2e1 + mrSges(3,2) * t527) * t358 + t527 * t709 + t275 * t51 + (t688 - t75 / 0.2e1) * Ifges(3,4)) * t676 + (-t48 * t193 + t574 * t33 + (-0.8e1 * (t27 * t530 + t36 * t697) * t358 + t210 * t742 + t51 * t738 + (0.2e1 * t75 + t551) * t334) * t350 * t114) * t359 + Ifges(3,4) * t114 * t551 + (t33 * t266 * t350 + (mrSges(3,2) * t728 + t523) * t739) * t358 + (-Ifges(2,6) * t48 + t33 * t624) * t350 + (Ifges(3,4) * t75 + t752 * t51) * t739 + (-t180 + (-mrSges(1,1) - t423) * t168) * t360 + t351 * (g(3) * t423 - t168 * t195 + t363);
t2 = -t92 * t17 - t122 * t12 - t125 * t14 + 0.8e1 * ((-t524 / 0.2e1 + mrSges(3,2) * t528) * t355 + t528 * t710 + t275 * t50 + (t689 - t74 / 0.2e1) * Ifges(3,4)) * t678 + (-t47 * t193 + t575 * t32 + (-0.8e1 * (t26 * t530 + t35 * t698) * t355 + t208 * t743 + t50 * t738 + (0.2e1 * t74 + t552) * t334) * t347 * t113) * t356 + Ifges(3,4) * t113 * t552 + (t32 * t266 * t347 + (mrSges(3,2) * t729 + t524) * t740) * t355 + (-Ifges(2,6) * t47 + t32 * t625) * t347 + (Ifges(3,4) * t74 + t753 * t50) * t740 + (-t180 + (-mrSges(1,1) - t424) * t167) * t357 + t348 * (g(3) * t424 - t167 * t195 + t363);
t1 = -t91 * t16 - t121 * t10 - t124 * t13 + 0.8e1 * ((-t525 / 0.2e1 + mrSges(3,2) * t529) * t352 + t529 * t711 + t275 * t49 + (t690 - t73 / 0.2e1) * Ifges(3,4)) * t680 + (-t46 * t193 + t576 * t31 + (-0.8e1 * (t25 * t530 + t34 * t699) * t352 + t206 * t744 + t49 * t738 + (0.2e1 * t73 + t553) * t334) * t344 * t112) * t353 + Ifges(3,4) * t112 * t553 + (t31 * t266 * t344 + (mrSges(3,2) * t730 + t525) * t741) * t352 + (-Ifges(2,6) * t46 + t31 * t626) * t344 + (Ifges(3,4) * t73 + t49 * t754) * t741 + (-t180 + (-mrSges(1,1) - t425) * t166) * t354 + t345 * (g(3) * t425 - t166 * t195 + t363);
t7 = [-g(1) * m(4) + t1 * t505 + t116 * t455 + t118 * t454 + t120 * t453 + t2 * t503 + t3 * t501 - t4 * t451 - t6 * t447 - t5 * t449 + (-t42 * t502 + (-t60 * t692 + t66 * t665) * t618 - t41 * t504 + (-t59 * t694 + t65 * t669) * t620 - t40 * t506 + (-t58 * t696 + t64 * t673) * t622) * t341 + (-t40 * t644 - t41 * t642 - t42 * t640 + (t482 * t60 + t483 * t59 + t484 * t58 + t510 * t66 + t511 * t65 + t512 * t64) * t406) * t340 + (t42 * t501 + (-t60 * t691 + t66 * t663) * t618 + t41 * t503 + (-t59 * t693 + t65 * t667) * t620 + t40 * t505 + (-t58 * t695 + t64 * t671) * t622 + m(4)) * t342; -g(2) * m(4) - t1 * t506 + t115 * t455 + t117 * t454 + t119 * t453 - t2 * t504 - t3 * t502 - t4 * t452 - t6 * t448 - t5 * t450 + (t39 * t501 + (-t57 * t691 + t63 * t663) * t618 + t38 * t503 + (-t56 * t693 + t62 * t667) * t620 + t37 * t505 + (-t55 * t695 + t61 * t671) * t622) * t342 + (-t37 * t644 - t38 * t642 - t39 * t640 + (t482 * t57 + t483 * t56 + t484 * t55 + t510 * t63 + t511 * t62 + t512 * t61) * t406) * t340 + (-t39 * t502 + (-t57 * t692 + t63 * t665) * t618 - t38 * t504 + (-t56 * t694 + t62 * t669) * t620 - t37 * t506 + (-t55 * t696 + t61 * t673) * t622 + m(4)) * t341; -g(3) * m(4) - t1 * t644 + t136 * t455 + t137 * t454 + t138 * t453 - t2 * t642 - t3 * t640 + t4 * t437 + t6 * t435 + t5 * t436 + (t54 * t501 + (t663 * t72 - t69 * t691) * t618 + t53 * t503 + (t667 * t71 - t68 * t693) * t620 + t52 * t505 + (-t67 * t695 + t671 * t70) * t622) * t342 + (-t54 * t502 + (t665 * t72 - t69 * t692) * t618 - t53 * t504 + (t669 * t71 - t68 * t694) * t620 - t52 * t506 + (-t67 * t696 + t673 * t70) * t622) * t341 + (-t52 * t644 - t53 * t642 - t54 * t640 + m(4) + (t482 * t69 + t483 * t68 + t484 * t67 + t510 * t72 + t511 * t71 + t512 * t70) * t406) * t340;];
tauX  = t7;
