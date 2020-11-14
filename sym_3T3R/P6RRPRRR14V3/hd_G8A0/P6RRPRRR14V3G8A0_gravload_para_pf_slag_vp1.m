% Calculate Gravitation load for parallel robot
% P6RRPRRR14V3G8A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [6x1]
%   Generalized platform coordinates
% qJ [3x6]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [6x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [6x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [6x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-12 23:36
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P6RRPRRR14V3G8A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,6),zeros(3,1),zeros(6,3),zeros(6,3),zeros(1,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [6 1]), ...
  'P6RRPRRR14V3G8A0_gravload_para_pf_slag_vp1: xP has to be [6x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 6]), ...
  'P6RRPRRR14V3G8A0_gravload_para_pf_slag_vp1: qJ has to be [3x6] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'P6RRPRRR14V3G8A0_gravload_para_pf_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P6RRPRRR14V3G8A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P6RRPRRR14V3G8A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P6RRPRRR14V3G8A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [6 3]), ...
  'P6RRPRRR14V3G8A0_gravload_para_pf_slag_vp1: legFrame has to be [6x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6RRPRRR14V3G8A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [6x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 1
% StartTime: 2020-03-12 23:32:40
% EndTime: 2020-03-12 23:32:41
% DurationCPUTime: 0.54s
% Computational Cost: add. (4451->481), mult. (9987->1076), div. (162->12), fcn. (9905->66), ass. (0->360)
unknown=NaN(6,1);
t1 = 0.1e1 / qJ(3,1);
t2 = sin(qJ(1,1));
t3 = cos(legFrame(1,3));
t5 = cos(qJ(1,1));
t6 = sin(legFrame(1,3));
t8 = t3 * t2 + t6 * t5;
t10 = cos(legFrame(1,2));
t11 = sin(qJ(2,1));
t12 = 0.1e1 / t11;
t16 = sin(legFrame(1,1));
t17 = sin(legFrame(1,2));
t18 = t17 * t16;
t20 = cos(legFrame(1,1));
t24 = t17 * t20;
t29 = -g(1) * t6 * t10 + g(2) * (-t6 * t18 + t3 * t20) + g(3) * (t3 * t16 + t6 * t24);
t30 = rSges(3,3) + qJ(3,1);
t32 = m(2) * rSges(2,2);
t33 = m(3) * t30 - t32;
t38 = rSges(2,1) * m(2) + m(3) * rSges(3,1);
t39 = t38 * t29;
t40 = cos(qJ(2,1));
t45 = m(1) * rSges(1,2) - m(2) * rSges(2,3) - m(3) * rSges(3,2);
t56 = g(1) * t3 * t10 + g(2) * (t3 * t18 + t6 * t20) + g(3) * (t6 * t16 - t3 * t24);
t64 = t38 * t56;
t71 = t5 * (-rSges(1,1) * t29 * m(1) - t11 * t33 * t29 - t40 * t39 + t56 * t45) + t2 * (rSges(1,1) * t56 * m(1) + t11 * t56 * t33 + t45 * t29 + t40 * t64);
t76 = -t2 * t6 + t5 * t3;
t80 = t10 * t76 * t40 + t17 * t11;
t83 = -m(3) * t30 + t32;
t89 = t10 * t16;
t91 = t10 * t20;
t93 = g(1) * t17 - g(2) * t89 + g(3) * t91;
t102 = t40 * (t2 * t83 * t29 + t5 * t83 * t56 - t38 * t93) + (t2 * t39 + t5 * t64 + t83 * t93) * t11;
t107 = t10 * t76 * t11 - t40 * t17;
t114 = t11 * t2 * t29 + t11 * t5 * t56 - t40 * t93;
t116 = 0.1e1 / qJ(3,2);
t117 = sin(qJ(1,2));
t118 = cos(legFrame(2,3));
t120 = cos(qJ(1,2));
t121 = sin(legFrame(2,3));
t123 = t118 * t117 + t121 * t120;
t125 = cos(legFrame(2,2));
t126 = sin(qJ(2,2));
t127 = 0.1e1 / t126;
t131 = sin(legFrame(2,1));
t132 = sin(legFrame(2,2));
t133 = t132 * t131;
t135 = cos(legFrame(2,1));
t139 = t132 * t135;
t144 = -g(1) * t121 * t125 + g(2) * (t118 * t135 - t121 * t133) + g(3) * (t118 * t131 + t121 * t139);
t145 = rSges(3,3) + qJ(3,2);
t147 = m(3) * t145 - t32;
t150 = t38 * t144;
t151 = cos(qJ(2,2));
t163 = g(1) * t118 * t125 + g(2) * (t118 * t133 + t121 * t135) + g(3) * (-t118 * t139 + t121 * t131);
t171 = t38 * t163;
t178 = t120 * (-rSges(1,1) * t144 * m(1) - t126 * t147 * t144 - t151 * t150 + t163 * t45) + t117 * (rSges(1,1) * t163 * m(1) + t126 * t163 * t147 + t45 * t144 + t151 * t171);
t183 = -t117 * t121 + t120 * t118;
t187 = t125 * t183 * t151 + t132 * t126;
t190 = -m(3) * t145 + t32;
t196 = t125 * t131;
t198 = t125 * t135;
t200 = g(1) * t132 - g(2) * t196 + g(3) * t198;
t209 = t151 * (t117 * t190 * t144 + t120 * t190 * t163 - t38 * t200) + (t117 * t150 + t120 * t171 + t190 * t200) * t126;
t214 = t125 * t183 * t126 - t151 * t132;
t221 = t126 * t117 * t144 + t126 * t120 * t163 - t151 * t200;
t223 = 0.1e1 / qJ(3,3);
t224 = sin(qJ(1,3));
t225 = cos(legFrame(3,3));
t227 = cos(qJ(1,3));
t228 = sin(legFrame(3,3));
t230 = t225 * t224 + t228 * t227;
t232 = cos(legFrame(3,2));
t233 = sin(qJ(2,3));
t234 = 0.1e1 / t233;
t238 = sin(legFrame(3,1));
t239 = sin(legFrame(3,2));
t240 = t239 * t238;
t242 = cos(legFrame(3,1));
t246 = t239 * t242;
t251 = -g(1) * t228 * t232 + g(2) * (t225 * t242 - t228 * t240) + g(3) * (t225 * t238 + t228 * t246);
t252 = rSges(3,3) + qJ(3,3);
t254 = m(3) * t252 - t32;
t257 = t38 * t251;
t258 = cos(qJ(2,3));
t270 = g(1) * t225 * t232 + g(2) * (t225 * t240 + t228 * t242) + g(3) * (-t225 * t246 + t228 * t238);
t278 = t38 * t270;
t285 = t227 * (-rSges(1,1) * t251 * m(1) - t233 * t254 * t251 - t258 * t257 + t270 * t45) + t224 * (rSges(1,1) * t270 * m(1) + t233 * t270 * t254 + t45 * t251 + t258 * t278);
t290 = -t224 * t228 + t227 * t225;
t294 = t232 * t290 * t258 + t239 * t233;
t297 = -m(3) * t252 + t32;
t303 = t232 * t238;
t305 = t232 * t242;
t307 = g(1) * t239 - g(2) * t303 + g(3) * t305;
t316 = t258 * (t224 * t297 * t251 + t227 * t297 * t270 - t38 * t307) + (t224 * t257 + t227 * t278 + t297 * t307) * t233;
t321 = t232 * t290 * t233 - t258 * t239;
t328 = t233 * t224 * t251 + t233 * t227 * t270 - t258 * t307;
t330 = 0.1e1 / qJ(3,4);
t331 = sin(qJ(1,4));
t332 = cos(legFrame(4,3));
t334 = cos(qJ(1,4));
t335 = sin(legFrame(4,3));
t337 = t332 * t331 + t335 * t334;
t339 = cos(legFrame(4,2));
t340 = sin(qJ(2,4));
t341 = 0.1e1 / t340;
t345 = sin(legFrame(4,1));
t346 = sin(legFrame(4,2));
t347 = t346 * t345;
t349 = cos(legFrame(4,1));
t353 = t346 * t349;
t358 = -g(1) * t335 * t339 + g(2) * (t332 * t349 - t335 * t347) + g(3) * (t332 * t345 + t335 * t353);
t359 = rSges(3,3) + qJ(3,4);
t361 = m(3) * t359 - t32;
t364 = t38 * t358;
t365 = cos(qJ(2,4));
t377 = g(1) * t332 * t339 + g(2) * (t332 * t347 + t335 * t349) + g(3) * (-t332 * t353 + t335 * t345);
t385 = t38 * t377;
t392 = t334 * (-rSges(1,1) * t358 * m(1) - t340 * t361 * t358 - t365 * t364 + t377 * t45) + t331 * (rSges(1,1) * t377 * m(1) + t340 * t377 * t361 + t45 * t358 + t365 * t385);
t397 = -t331 * t335 + t334 * t332;
t401 = t339 * t397 * t365 + t346 * t340;
t404 = -m(3) * t359 + t32;
t410 = t339 * t345;
t412 = t339 * t349;
t414 = g(1) * t346 - g(2) * t410 + g(3) * t412;
t423 = t365 * (t331 * t404 * t358 + t334 * t404 * t377 - t38 * t414) + (t331 * t364 + t334 * t385 + t404 * t414) * t340;
t428 = t339 * t397 * t340 - t365 * t346;
t435 = t340 * t331 * t358 + t340 * t334 * t377 - t365 * t414;
t437 = 0.1e1 / qJ(3,5);
t438 = sin(qJ(1,5));
t439 = cos(legFrame(5,3));
t441 = cos(qJ(1,5));
t442 = sin(legFrame(5,3));
t444 = t439 * t438 + t442 * t441;
t446 = cos(legFrame(5,2));
t447 = sin(qJ(2,5));
t448 = 0.1e1 / t447;
t452 = sin(legFrame(5,1));
t453 = sin(legFrame(5,2));
t454 = t453 * t452;
t456 = cos(legFrame(5,1));
t460 = t453 * t456;
t465 = -g(1) * t442 * t446 + g(2) * (t439 * t456 - t442 * t454) + g(3) * (t439 * t452 + t442 * t460);
t466 = rSges(3,3) + qJ(3,5);
t468 = m(3) * t466 - t32;
t471 = t38 * t465;
t472 = cos(qJ(2,5));
t484 = g(1) * t439 * t446 + g(2) * (t439 * t454 + t442 * t456) + g(3) * (-t439 * t460 + t442 * t452);
t492 = t38 * t484;
t499 = t441 * (-rSges(1,1) * t465 * m(1) - t447 * t468 * t465 + t484 * t45 - t472 * t471) + t438 * (rSges(1,1) * t484 * m(1) + t447 * t484 * t468 + t45 * t465 + t472 * t492);
t504 = -t438 * t442 + t441 * t439;
t508 = t446 * t504 * t472 + t453 * t447;
t511 = -m(3) * t466 + t32;
t517 = t446 * t452;
t519 = t446 * t456;
t521 = g(1) * t453 - g(2) * t517 + g(3) * t519;
t530 = t472 * (t438 * t511 * t465 + t441 * t511 * t484 - t38 * t521) + (t438 * t471 + t441 * t492 + t511 * t521) * t447;
t535 = t446 * t504 * t447 - t472 * t453;
t542 = t447 * t438 * t465 + t447 * t441 * t484 - t472 * t521;
t544 = 0.1e1 / qJ(3,6);
t545 = sin(qJ(1,6));
t546 = cos(legFrame(6,3));
t548 = cos(qJ(1,6));
t549 = sin(legFrame(6,3));
t551 = t546 * t545 + t549 * t548;
t553 = cos(legFrame(6,2));
t554 = sin(qJ(2,6));
t555 = 0.1e1 / t554;
t559 = sin(legFrame(6,1));
t560 = sin(legFrame(6,2));
t561 = t560 * t559;
t563 = cos(legFrame(6,1));
t567 = t560 * t563;
t572 = -g(1) * t549 * t553 + g(2) * (t546 * t563 - t549 * t561) + g(3) * (t546 * t559 + t549 * t567);
t573 = rSges(3,3) + qJ(3,6);
t575 = m(3) * t573 - t32;
t578 = t38 * t572;
t579 = cos(qJ(2,6));
t591 = g(1) * t546 * t553 + g(2) * (t546 * t561 + t549 * t563) + g(3) * (-t546 * t567 + t549 * t559);
t599 = t38 * t591;
t606 = t548 * (-rSges(1,1) * t572 * m(1) - t554 * t575 * t572 + t591 * t45 - t579 * t578) + t545 * (rSges(1,1) * t591 * m(1) + t554 * t591 * t575 + t45 * t572 + t579 * t599);
t611 = -t545 * t549 + t548 * t546;
t615 = t553 * t611 * t579 + t560 * t554;
t618 = -m(3) * t573 + t32;
t624 = t553 * t559;
t626 = t553 * t563;
t628 = g(1) * t560 - g(2) * t624 + g(3) * t626;
t637 = t579 * (t545 * t618 * t572 + t548 * t618 * t591 - t38 * t628) + (t545 * t578 + t548 * t599 + t618 * t628) * t554;
t642 = t553 * t611 * t554 - t579 * t560;
t649 = t554 * t545 * t572 + t554 * t548 * t591 - t579 * t628;
t652 = -t71 * t12 * t10 * t8 * t1 + t102 * t1 * t80 - t114 * m(3) * t107 - t178 * t127 * t125 * t123 * t116 + t209 * t116 * t187 - t221 * m(3) * t214 - t285 * t234 * t232 * t230 * t223 + t316 * t223 * t294 - t328 * m(3) * t321 - t392 * t341 * t339 * t337 * t330 + t423 * t330 * t401 - t435 * m(3) * t428 - t499 * t448 * t446 * t444 * t437 + t530 * t437 * t508 - t542 * m(3) * t535 - t606 * t555 * t553 * t551 * t544 + t637 * t544 * t615 - t649 * m(3) * t642 - m(4) * g(1);
t655 = -t8 * t18 + t20 * t76;
t657 = t71 * t12;
t659 = t76 * t17;
t662 = t16 * t659 + t8 * t20;
t666 = -t10 * t11 * t16 + t40 * t662;
t671 = t11 * t662 + t40 * t89;
t676 = -t123 * t133 + t135 * t183;
t678 = t178 * t127;
t680 = t183 * t132;
t683 = t123 * t135 + t131 * t680;
t687 = -t125 * t126 * t131 + t151 * t683;
t692 = t126 * t683 + t151 * t196;
t697 = -t230 * t240 + t242 * t290;
t699 = t285 * t234;
t701 = t290 * t239;
t704 = t230 * t242 + t238 * t701;
t708 = -t232 * t233 * t238 + t258 * t704;
t713 = t233 * t704 + t258 * t303;
t718 = -t337 * t347 + t349 * t397;
t720 = t392 * t341;
t722 = t397 * t346;
t725 = t337 * t349 + t345 * t722;
t729 = -t339 * t340 * t345 + t365 * t725;
t734 = t340 * t725 + t365 * t410;
t739 = -t444 * t454 + t456 * t504;
t741 = t499 * t448;
t743 = t504 * t453;
t746 = t444 * t456 + t452 * t743;
t750 = -t446 * t447 * t452 + t472 * t746;
t755 = t447 * t746 + t472 * t517;
t760 = -t551 * t561 + t563 * t611;
t762 = t606 * t555;
t764 = t611 * t560;
t767 = t551 * t563 + t559 * t764;
t771 = -t553 * t554 * t559 + t579 * t767;
t776 = t554 * t767 + t579 * t624;
t780 = t657 * t1 * t655 + t102 * t1 * t666 - t114 * m(3) * t671 + t678 * t116 * t676 + t209 * t116 * t687 - t221 * m(3) * t692 + t699 * t223 * t697 + t316 * t223 * t708 - t328 * m(3) * t713 + t720 * t330 * t718 + t423 * t330 * t729 - t435 * m(3) * t734 + t741 * t437 * t739 + t530 * t437 * t750 - t542 * m(3) * t755 + t762 * t544 * t760 + t637 * t544 * t771 - t649 * m(3) * t776 - m(4) * g(2);
t784 = t20 * t8 * t17 + t76 * t16;
t789 = t8 * t16 - t20 * t659;
t793 = t10 * t11 * t20 + t40 * t789;
t798 = t11 * t789 - t40 * t91;
t804 = t135 * t123 * t132 + t183 * t131;
t809 = t123 * t131 - t135 * t680;
t813 = t125 * t126 * t135 + t151 * t809;
t818 = t126 * t809 - t151 * t198;
t824 = t242 * t230 * t239 + t290 * t238;
t829 = t230 * t238 - t242 * t701;
t833 = t232 * t233 * t242 + t258 * t829;
t838 = t233 * t829 - t258 * t305;
t844 = t349 * t337 * t346 + t397 * t345;
t849 = t337 * t345 - t349 * t722;
t853 = t339 * t340 * t349 + t365 * t849;
t858 = t340 * t849 - t365 * t412;
t864 = t456 * t444 * t453 + t504 * t452;
t869 = t444 * t452 - t456 * t743;
t873 = t446 * t447 * t456 + t472 * t869;
t878 = t447 * t869 - t472 * t519;
t884 = t563 * t551 * t560 + t611 * t559;
t889 = t551 * t559 - t563 * t764;
t893 = t553 * t554 * t563 + t579 * t889;
t898 = t554 * t889 - t579 * t626;
t902 = t657 * t1 * t784 + t102 * t1 * t793 - t114 * m(3) * t798 + t678 * t116 * t804 + t209 * t116 * t813 - t221 * m(3) * t818 + t699 * t223 * t824 + t316 * t223 * t833 - t328 * m(3) * t838 + t720 * t330 * t844 + t423 * t330 * t853 - t435 * m(3) * t858 + t741 * t437 * t864 + t530 * t437 * t873 - t542 * m(3) * t878 + t762 * t544 * t884 + t637 * t544 * t893 - t649 * m(3) * t898 - m(4) * g(3);
t903 = sin(xP(6));
t905 = cos(xP(6));
t907 = -koppelP(1,2) * t903 + t905 * koppelP(1,1);
t908 = sin(xP(5));
t910 = cos(xP(5));
t911 = koppelP(1,3) * t910;
t912 = t908 * t907 - t911;
t913 = cos(xP(4));
t915 = sin(xP(4));
t918 = koppelP(1,1) * t903 + t905 * koppelP(1,2);
t919 = t918 * t915;
t920 = t913 * t912 - t919;
t922 = t12 * t1;
t925 = t918 * t913;
t926 = t915 * t912 + t925;
t944 = -koppelP(2,2) * t903 + t905 * koppelP(2,1);
t946 = koppelP(2,3) * t910;
t947 = t908 * t944 - t946;
t951 = koppelP(2,1) * t903 + t905 * koppelP(2,2);
t952 = t951 * t915;
t953 = t913 * t947 - t952;
t955 = t127 * t116;
t958 = t951 * t913;
t959 = t915 * t947 + t958;
t977 = -koppelP(3,2) * t903 + t905 * koppelP(3,1);
t979 = koppelP(3,3) * t910;
t980 = t908 * t977 - t979;
t984 = koppelP(3,1) * t903 + t905 * koppelP(3,2);
t985 = t984 * t915;
t986 = t913 * t980 - t985;
t988 = t234 * t223;
t991 = t984 * t913;
t992 = t915 * t980 + t991;
t1010 = -koppelP(4,2) * t903 + t905 * koppelP(4,1);
t1012 = koppelP(4,3) * t910;
t1013 = t908 * t1010 - t1012;
t1017 = koppelP(4,1) * t903 + t905 * koppelP(4,2);
t1018 = t1017 * t915;
t1019 = t913 * t1013 - t1018;
t1021 = t341 * t330;
t1024 = t1017 * t913;
t1025 = t915 * t1013 + t1024;
t1043 = -koppelP(5,2) * t903 + t905 * koppelP(5,1);
t1045 = koppelP(5,3) * t910;
t1046 = t908 * t1043 - t1045;
t1050 = koppelP(5,1) * t903 + t905 * koppelP(5,2);
t1051 = t1050 * t915;
t1052 = t913 * t1046 - t1051;
t1054 = t448 * t437;
t1057 = t1050 * t913;
t1058 = t915 * t1046 + t1057;
t1076 = -koppelP(6,2) * t903 + t905 * koppelP(6,1);
t1078 = koppelP(6,3) * t910;
t1079 = t908 * t1076 - t1078;
t1083 = koppelP(6,1) * t903 + t905 * koppelP(6,2);
t1084 = t1083 * t915;
t1085 = t913 * t1079 - t1084;
t1087 = t555 * t544;
t1090 = t1083 * t913;
t1091 = t915 * t1079 + t1090;
t1107 = t908 * g(2);
t1109 = g(3) * rSges(4,2);
t1112 = rSges(4,2) * t903;
t1120 = g(3) * t908;
t1127 = t903 * g(2);
t1135 = t71 * (t922 * t655 * t920 + t922 * t784 * t926) + t102 * (t1 * t666 * t920 + t1 * t793 * t926) - t114 * m(3) * (t671 * t920 + t798 * t926) + t178 * (t955 * t676 * t953 + t955 * t804 * t959) + t209 * (t116 * t687 * t953 + t116 * t813 * t959) - t221 * m(3) * (t692 * t953 + t818 * t959) + t285 * (t988 * t697 * t986 + t988 * t824 * t992) + t316 * (t223 * t708 * t986 + t223 * t833 * t992) - t328 * m(3) * (t713 * t986 + t838 * t992) + t392 * (t1021 * t718 * t1019 + t1021 * t844 * t1025) + t423 * (t330 * t729 * t1019 + t330 * t853 * t1025) - t435 * m(3) * (t734 * t1019 + t858 * t1025) + t499 * (t1054 * t739 * t1052 + t1054 * t864 * t1058) + t530 * (t437 * t750 * t1052 + t437 * t873 * t1058) - t542 * m(3) * (t755 * t1052 + t878 * t1058) + t606 * (t1087 * t760 * t1085 + t1087 * t884 * t1091) + t637 * (t544 * t771 * t1085 + t544 * t893 * t1091) - t649 * m(3) * (t776 * t1085 + t898 * t1091) + m(4) * (t913 * (t905 * (-rSges(4,1) * t1107 - t1109) + t1112 * t1107 + rSges(4,3) * t910 * g(2) - t903 * rSges(4,1) * g(3)) + t915 * (t905 * (g(2) * rSges(4,2) - rSges(4,1) * t1120) + t1109 * t903 * t908 + rSges(4,1) * t1127 + rSges(4,3) * g(3) * t910));
t1137 = -t908 * t907 + t911;
t1139 = t913 * t1137 + t919;
t1142 = t12 * t10 * t8;
t1144 = t905 * t910;
t1146 = t903 * t910;
t1149 = -koppelP(1,1) * t1144 + koppelP(1,2) * t1146 - koppelP(1,3) * t908;
t1166 = -t908 * t944 + t946;
t1168 = t913 * t1166 + t952;
t1171 = t127 * t125 * t123;
t1176 = -koppelP(2,1) * t1144 + koppelP(2,2) * t1146 - koppelP(2,3) * t908;
t1193 = -t908 * t977 + t979;
t1195 = t913 * t1193 + t985;
t1198 = t234 * t232 * t230;
t1203 = -koppelP(3,1) * t1144 + koppelP(3,2) * t1146 - koppelP(3,3) * t908;
t1220 = -t908 * t1010 + t1012;
t1222 = t913 * t1220 + t1018;
t1225 = t341 * t339 * t337;
t1230 = -koppelP(4,1) * t1144 + koppelP(4,2) * t1146 - koppelP(4,3) * t908;
t1247 = -t908 * t1043 + t1045;
t1249 = t913 * t1247 + t1051;
t1252 = t448 * t446 * t444;
t1257 = -koppelP(5,1) * t1144 + koppelP(5,2) * t1146 - koppelP(5,3) * t908;
t1274 = -t908 * t1076 + t1078;
t1276 = t913 * t1274 + t1084;
t1279 = t555 * t553 * t551;
t1284 = -koppelP(6,1) * t1144 + koppelP(6,2) * t1146 - koppelP(6,3) * t908;
t1302 = rSges(4,1) * t905 - t1112;
t1313 = g(1) * rSges(4,1);
t1318 = t71 * (-t1142 * t1 * t1139 + t922 * t784 * t1149) + t102 * (t1 * t80 * t1139 + t1 * t793 * t1149) - t114 * m(3) * (t107 * t1139 + t798 * t1149) + t178 * (-t1171 * t116 * t1168 + t955 * t804 * t1176) + t209 * (t116 * t187 * t1168 + t116 * t813 * t1176) - t221 * m(3) * (t214 * t1168 + t818 * t1176) + t285 * (-t1198 * t223 * t1195 + t988 * t824 * t1203) + t316 * (t223 * t294 * t1195 + t223 * t833 * t1203) - t328 * m(3) * (t321 * t1195 + t838 * t1203) + t392 * (t1021 * t844 * t1230 - t1225 * t330 * t1222) + t423 * (t330 * t401 * t1222 + t330 * t853 * t1230) - t435 * m(3) * (t428 * t1222 + t858 * t1230) + t499 * (t1054 * t864 * t1257 - t1252 * t437 * t1249) + t530 * (t437 * t508 * t1249 + t437 * t873 * t1257) - t542 * m(3) * (t535 * t1249 + t878 * t1257) + t606 * (t1087 * t884 * t1284 - t1279 * t544 * t1276) + t637 * (t544 * t615 * t1276 + t544 * t893 * t1284) - t649 * m(3) * (t642 * t1276 + t898 * t1284) - m(4) * (t913 * g(1) * (-t1302 * t908 + rSges(4,3) * t910) - t910 * t1302 * g(3) + g(1) * rSges(4,2) * t915 * t905 + t1313 * t903 * t915 - rSges(4,3) * t1120);
t1320 = t915 * t1137 - t925;
t1339 = t915 * t1166 - t958;
t1358 = t915 * t1193 - t991;
t1377 = t915 * t1220 - t1024;
t1396 = t915 * t1247 - t1057;
t1415 = t915 * t1274 - t1090;
t1455 = t71 * (-t1142 * t1 * t1320 - t922 * t655 * t1149) + t102 * (-t1 * t666 * t1149 + t1 * t80 * t1320) - t114 * m(3) * (t107 * t1320 - t671 * t1149) + t178 * (-t1171 * t116 * t1339 - t955 * t676 * t1176) + t209 * (-t116 * t687 * t1176 + t116 * t187 * t1339) - t221 * m(3) * (-t692 * t1176 + t214 * t1339) + t285 * (-t1198 * t223 * t1358 - t988 * t697 * t1203) + t316 * (-t223 * t708 * t1203 + t223 * t294 * t1358) - t328 * m(3) * (-t713 * t1203 + t321 * t1358) + t392 * (-t1021 * t718 * t1230 - t1225 * t330 * t1377) + t423 * (-t330 * t729 * t1230 + t330 * t401 * t1377) - t435 * m(3) * (-t734 * t1230 + t428 * t1377) + t499 * (-t1054 * t739 * t1257 - t1252 * t437 * t1396) + t530 * (-t437 * t750 * t1257 + t437 * t508 * t1396) - t542 * m(3) * (-t755 * t1257 + t535 * t1396) + t606 * (-t1087 * t760 * t1284 - t1279 * t544 * t1415) + t637 * (-t544 * t771 * t1284 + t544 * t615 * t1415) - t649 * m(3) * (-t776 * t1284 + t642 * t1415) + m(4) * (t910 * (-rSges(4,3) * t915 * g(1) - t905 * rSges(4,1) * g(2) + rSges(4,2) * t1127) + t905 * (t915 * rSges(4,1) * t908 + rSges(4,2) * t913) * g(1) - rSges(4,2) * g(1) * t903 * t915 * t908 + t1313 * t903 * t913 - rSges(4,3) * t1107);
unknown(1,1) = t652;
unknown(2,1) = t780;
unknown(3,1) = t902;
unknown(4,1) = t1135;
unknown(5,1) = t1318;
unknown(6,1) = t1455;
taugX  = unknown;
