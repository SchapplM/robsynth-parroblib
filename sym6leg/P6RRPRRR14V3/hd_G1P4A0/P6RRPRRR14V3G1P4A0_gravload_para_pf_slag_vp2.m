% Calculate Gravitation load for parallel robot
% P6RRPRRR14V3G1P4A0
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
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [6x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-12 23:28
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P6RRPRRR14V3G1P4A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,6),zeros(3,1),zeros(6,3),zeros(6,3),zeros(1,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [6 1]), ...
  'P6RRPRRR14V3G1P4A0_gravload_para_pf_slag_vp2: xP has to be [6x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 6]), ...
  'P6RRPRRR14V3G1P4A0_gravload_para_pf_slag_vp2: qJ has to be [3x6] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'P6RRPRRR14V3G1P4A0_gravload_para_pf_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P6RRPRRR14V3G1P4A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P6RRPRRR14V3G1P4A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P6RRPRRR14V3G1P4A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [6 3]), ...
  'P6RRPRRR14V3G1P4A0_gravload_para_pf_slag_vp2: legFrame has to be [6x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6RRPRRR14V3G1P4A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [6x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 1
% StartTime: 2020-03-12 23:27:16
% EndTime: 2020-03-12 23:27:16
% DurationCPUTime: 0.31s
% Computational Cost: add. (2141->355), mult. (3378->737), div. (126->12), fcn. (2783->42), ass. (0->273)
unknown=NaN(6,1);
t1 = sin(legFrame(1,3));
t2 = cos(qJ(1,1));
t4 = sin(qJ(1,1));
t5 = cos(legFrame(1,3));
t7 = -t2 * t1 - t5 * t4;
t8 = sin(qJ(2,1));
t9 = 0.1e1 / t8;
t11 = 0.1e1 / qJ(3,1);
t14 = -g(1) * t1 + g(2) * t5;
t16 = m(3) * qJ(3,1) - mrSges(2,2) + mrSges(3,3);
t17 = t16 * t14;
t19 = mrSges(2,1) + mrSges(3,1);
t20 = t14 * t19;
t21 = cos(qJ(2,1));
t23 = mrSges(1,2) - mrSges(3,2) - mrSges(2,3);
t26 = g(1) * t5 + g(2) * t1;
t31 = t16 * t26;
t33 = t19 * t26;
t39 = t2 * (-mrSges(1,1) * t14 - t8 * t17 - t21 * t20 + t26 * t23) + t4 * (mrSges(1,1) * t26 + t23 * t14 + t21 * t33 + t8 * t31);
t40 = t39 * t11;
t44 = -t4 * t1 + t5 * t2;
t48 = t19 * g(3);
t56 = t21 * (-t4 * t17 - t2 * t31 - t48) - t8 * (t16 * g(3) - t2 * t33 - t4 * t20);
t65 = t14 * t4 * t8 + t26 * t8 * t2 - g(3) * t21;
t66 = t65 * m(3);
t68 = sin(legFrame(2,3));
t69 = cos(qJ(1,2));
t71 = sin(qJ(1,2));
t72 = cos(legFrame(2,3));
t74 = -t69 * t68 - t72 * t71;
t75 = sin(qJ(2,2));
t76 = 0.1e1 / t75;
t78 = 0.1e1 / qJ(3,2);
t81 = -g(1) * t68 + g(2) * t72;
t83 = m(3) * qJ(3,2) - mrSges(2,2) + mrSges(3,3);
t84 = t83 * t81;
t86 = t19 * t81;
t87 = cos(qJ(2,2));
t91 = g(1) * t72 + g(2) * t68;
t96 = t83 * t91;
t98 = t19 * t91;
t104 = t69 * (-mrSges(1,1) * t81 + t91 * t23 - t75 * t84 - t87 * t86) + t71 * (mrSges(1,1) * t91 + t23 * t81 + t75 * t96 + t87 * t98);
t105 = t104 * t78;
t109 = -t71 * t68 + t72 * t69;
t120 = t87 * (-t69 * t96 - t71 * t84 - t48) - t75 * (t83 * g(3) - t69 * t98 - t71 * t86);
t129 = t91 * t75 * t69 + t81 * t71 * t75 - g(3) * t87;
t130 = t129 * m(3);
t132 = sin(legFrame(3,3));
t133 = cos(qJ(1,3));
t135 = sin(qJ(1,3));
t136 = cos(legFrame(3,3));
t138 = -t133 * t132 - t136 * t135;
t139 = sin(qJ(2,3));
t140 = 0.1e1 / t139;
t142 = 0.1e1 / qJ(3,3);
t145 = -g(1) * t132 + g(2) * t136;
t147 = m(3) * qJ(3,3) - mrSges(2,2) + mrSges(3,3);
t148 = t147 * t145;
t150 = t19 * t145;
t151 = cos(qJ(2,3));
t155 = g(1) * t136 + g(2) * t132;
t160 = t147 * t155;
t162 = t19 * t155;
t168 = t133 * (-mrSges(1,1) * t145 - t139 * t148 - t151 * t150 + t155 * t23) + t135 * (mrSges(1,1) * t155 + t139 * t160 + t23 * t145 + t151 * t162);
t169 = t168 * t142;
t173 = -t135 * t132 + t136 * t133;
t184 = t151 * (-t133 * t160 - t135 * t148 - t48) - t139 * (t147 * g(3) - t133 * t162 - t135 * t150);
t193 = t155 * t139 * t133 + t145 * t135 * t139 - g(3) * t151;
t194 = t193 * m(3);
t196 = sin(legFrame(4,3));
t197 = cos(qJ(1,4));
t199 = sin(qJ(1,4));
t200 = cos(legFrame(4,3));
t202 = -t197 * t196 - t200 * t199;
t203 = sin(qJ(2,4));
t204 = 0.1e1 / t203;
t206 = 0.1e1 / qJ(3,4);
t209 = -g(1) * t196 + g(2) * t200;
t211 = m(3) * qJ(3,4) - mrSges(2,2) + mrSges(3,3);
t212 = t211 * t209;
t214 = t19 * t209;
t215 = cos(qJ(2,4));
t219 = g(1) * t200 + g(2) * t196;
t224 = t211 * t219;
t226 = t19 * t219;
t232 = t197 * (-mrSges(1,1) * t209 - t203 * t212 - t215 * t214 + t219 * t23) + t199 * (mrSges(1,1) * t219 + t203 * t224 + t23 * t209 + t215 * t226);
t233 = t232 * t206;
t237 = -t199 * t196 + t200 * t197;
t248 = t215 * (-t197 * t224 - t199 * t212 - t48) - t203 * (t211 * g(3) - t197 * t226 - t199 * t214);
t257 = t219 * t203 * t197 + t209 * t199 * t203 - g(3) * t215;
t258 = t257 * m(3);
t260 = sin(legFrame(5,3));
t261 = cos(qJ(1,5));
t263 = sin(qJ(1,5));
t264 = cos(legFrame(5,3));
t266 = -t261 * t260 - t264 * t263;
t267 = sin(qJ(2,5));
t268 = 0.1e1 / t267;
t270 = 0.1e1 / qJ(3,5);
t273 = -g(1) * t260 + g(2) * t264;
t275 = m(3) * qJ(3,5) - mrSges(2,2) + mrSges(3,3);
t276 = t275 * t273;
t278 = t19 * t273;
t279 = cos(qJ(2,5));
t283 = g(1) * t264 + g(2) * t260;
t288 = t275 * t283;
t290 = t19 * t283;
t296 = t261 * (-mrSges(1,1) * t273 + t283 * t23 - t267 * t276 - t279 * t278) + t263 * (mrSges(1,1) * t283 + t23 * t273 + t267 * t288 + t279 * t290);
t297 = t296 * t270;
t301 = -t263 * t260 + t264 * t261;
t312 = t279 * (-t261 * t288 - t263 * t276 - t48) - t267 * (t275 * g(3) - t261 * t290 - t263 * t278);
t321 = t283 * t267 * t261 + t273 * t263 * t267 - g(3) * t279;
t322 = t321 * m(3);
t324 = sin(legFrame(6,3));
t325 = cos(qJ(1,6));
t327 = sin(qJ(1,6));
t328 = cos(legFrame(6,3));
t330 = -t325 * t324 - t328 * t327;
t331 = sin(qJ(2,6));
t332 = 0.1e1 / t331;
t334 = 0.1e1 / qJ(3,6);
t337 = -g(1) * t324 + g(2) * t328;
t339 = m(3) * qJ(3,6) - mrSges(2,2) + mrSges(3,3);
t340 = t339 * t337;
t342 = t19 * t337;
t343 = cos(qJ(2,6));
t347 = g(1) * t328 + g(2) * t324;
t352 = t339 * t347;
t354 = t19 * t347;
t360 = t325 * (-mrSges(1,1) * t337 + t347 * t23 - t331 * t340 - t343 * t342) + t327 * (mrSges(1,1) * t347 + t23 * t337 + t331 * t352 + t343 * t354);
t361 = t360 * t334;
t365 = -t327 * t324 + t328 * t325;
t376 = t343 * (-t325 * t352 - t327 * t340 - t48) - t331 * (t339 * g(3) - t325 * t354 - t327 * t342);
t385 = t347 * t331 * t325 + t337 * t327 * t331 - g(3) * t343;
t386 = t385 * m(3);
t389 = t40 * t9 * t7 + t56 * t11 * t21 * t44 - t66 * t8 * t44 + t105 * t76 * t74 + t120 * t78 * t87 * t109 - t130 * t75 * t109 + t169 * t140 * t138 + t184 * t142 * t151 * t173 - t194 * t139 * t173 + t233 * t204 * t202 + t248 * t206 * t215 * t237 - t258 * t203 * t237 + t297 * t268 * t266 + t312 * t270 * t279 * t301 - t322 * t267 * t301 + t361 * t332 * t330 + t376 * t334 * t343 * t365 - t386 * t331 * t365 - g(1) * m(4);
t433 = t40 * t9 * t44 - t56 * t21 * t11 * t7 + t66 * t8 * t7 + t105 * t76 * t109 - t120 * t87 * t74 * t78 + t130 * t74 * t75 + t169 * t140 * t173 - t184 * t151 * t138 * t142 + t194 * t138 * t139 + t233 * t204 * t237 - t248 * t215 * t202 * t206 + t258 * t202 * t203 + t297 * t268 * t301 - t312 * t279 * t266 * t270 + t322 * t266 * t267 + t361 * t332 * t365 - t376 * t343 * t330 * t334 + t386 * t330 * t331 - g(2) * m(4);
t459 = t56 * t8 * t11 + t65 * m(3) * t21 + t120 * t75 * t78 + t129 * m(3) * t87 + t184 * t139 * t142 + t193 * m(3) * t151 + t248 * t203 * t206 + t257 * m(3) * t215 + t312 * t267 * t270 + t321 * m(3) * t279 + t376 * t331 * t334 + t385 * m(3) * t343 - g(3) * m(4);
t460 = sin(xP(6));
t462 = cos(xP(6));
t464 = -koppelP(1,2) * t460 + t462 * koppelP(1,1);
t465 = sin(xP(5));
t467 = cos(xP(5));
t468 = koppelP(1,3) * t467;
t469 = t465 * t464 - t468;
t470 = cos(xP(4));
t472 = sin(xP(4));
t475 = koppelP(1,1) * t460 + koppelP(1,2) * t462;
t476 = t475 * t472;
t477 = t470 * t469 - t476;
t479 = t11 * t9;
t480 = t39 * t479;
t483 = -t21 * t7;
t486 = t475 * t470;
t487 = t472 * t469 + t486;
t500 = -koppelP(2,2) * t460 + t462 * koppelP(2,1);
t502 = koppelP(2,3) * t467;
t503 = t465 * t500 - t502;
t507 = koppelP(2,1) * t460 + koppelP(2,2) * t462;
t508 = t507 * t472;
t509 = t470 * t503 - t508;
t511 = t78 * t76;
t512 = t104 * t511;
t515 = -t87 * t74;
t518 = t507 * t470;
t519 = t472 * t503 + t518;
t532 = -koppelP(3,2) * t460 + t462 * koppelP(3,1);
t534 = koppelP(3,3) * t467;
t535 = t465 * t532 - t534;
t539 = koppelP(3,1) * t460 + koppelP(3,2) * t462;
t540 = t539 * t472;
t541 = t470 * t535 - t540;
t543 = t142 * t140;
t544 = t168 * t543;
t547 = -t151 * t138;
t550 = t539 * t470;
t551 = t472 * t535 + t550;
t564 = -koppelP(4,2) * t460 + t462 * koppelP(4,1);
t566 = koppelP(4,3) * t467;
t567 = t465 * t564 - t566;
t571 = koppelP(4,1) * t460 + koppelP(4,2) * t462;
t572 = t571 * t472;
t573 = t470 * t567 - t572;
t575 = t206 * t204;
t576 = t232 * t575;
t579 = -t215 * t202;
t582 = t571 * t470;
t583 = t472 * t567 + t582;
t596 = -koppelP(5,2) * t460 + t462 * koppelP(5,1);
t598 = koppelP(5,3) * t467;
t599 = t465 * t596 - t598;
t603 = koppelP(5,1) * t460 + koppelP(5,2) * t462;
t604 = t603 * t472;
t605 = t470 * t599 - t604;
t607 = t270 * t268;
t608 = t296 * t607;
t611 = -t279 * t266;
t614 = t603 * t470;
t615 = t472 * t599 + t614;
t628 = -koppelP(6,2) * t460 + t462 * koppelP(6,1);
t630 = koppelP(6,3) * t467;
t631 = t465 * t628 - t630;
t635 = koppelP(6,1) * t460 + t462 * koppelP(6,2);
t636 = t635 * t472;
t637 = t470 * t631 - t636;
t639 = t334 * t332;
t640 = t360 * t639;
t643 = -t343 * t330;
t646 = t635 * t470;
t647 = t472 * t631 + t646;
t658 = t465 * g(2);
t660 = g(3) * mrSges(4,2);
t663 = mrSges(4,2) * t460;
t671 = g(3) * t465;
t678 = t460 * g(2);
t684 = t480 * t44 * t477 + t56 * (t483 * t11 * t477 + t8 * t11 * t487) - t65 * m(3) * (-t7 * t8 * t477 - t21 * t487) + t512 * t109 * t509 + t120 * (t515 * t78 * t509 + t75 * t78 * t519) - t129 * m(3) * (-t74 * t75 * t509 - t87 * t519) + t544 * t173 * t541 + t184 * (t139 * t142 * t551 + t547 * t142 * t541) - t193 * m(3) * (-t138 * t139 * t541 - t151 * t551) + t576 * t237 * t573 + t248 * (t203 * t206 * t583 + t579 * t206 * t573) - t257 * m(3) * (-t202 * t203 * t573 - t215 * t583) + t608 * t301 * t605 + t312 * (t267 * t270 * t615 + t611 * t270 * t605) - t321 * m(3) * (-t266 * t267 * t605 - t279 * t615) + t640 * t365 * t637 + t376 * (t331 * t334 * t647 + t643 * t334 * t637) - t385 * m(3) * (-t330 * t331 * t637 - t343 * t647) - t470 * (t462 * (mrSges(4,1) * t658 + t660) - t663 * t658 - mrSges(4,3) * t467 * g(2) + t460 * mrSges(4,1) * g(3)) + t472 * (t462 * (g(2) * mrSges(4,2) - mrSges(4,1) * t671) + t660 * t460 * t465 + mrSges(4,1) * t678 + mrSges(4,3) * g(3) * t467);
t686 = -t465 * t464 + t468;
t688 = t470 * t686 + t476;
t691 = t44 * t688;
t692 = t11 * t21;
t694 = t462 * t467;
t696 = t460 * t467;
t699 = -koppelP(1,3) * t465 - koppelP(1,1) * t694 + koppelP(1,2) * t696;
t710 = -t465 * t500 + t502;
t712 = t470 * t710 + t508;
t715 = t109 * t712;
t716 = t78 * t87;
t721 = -koppelP(2,3) * t465 - koppelP(2,1) * t694 + koppelP(2,2) * t696;
t732 = -t465 * t532 + t534;
t734 = t470 * t732 + t540;
t737 = t173 * t734;
t738 = t142 * t151;
t743 = -koppelP(3,3) * t465 - koppelP(3,1) * t694 + koppelP(3,2) * t696;
t754 = -t465 * t564 + t566;
t756 = t470 * t754 + t572;
t759 = t237 * t756;
t760 = t206 * t215;
t765 = -koppelP(4,3) * t465 - koppelP(4,1) * t694 + koppelP(4,2) * t696;
t770 = t480 * t7 * t688 + t56 * (t8 * t11 * t699 + t692 * t691) - t65 * m(3) * (-t21 * t699 + t8 * t691) + t512 * t74 * t712 + t120 * (t75 * t78 * t721 + t716 * t715) - t129 * m(3) * (t75 * t715 - t87 * t721) + t544 * t138 * t734 + t184 * (t139 * t142 * t743 + t738 * t737) - t193 * m(3) * (t139 * t737 - t151 * t743) + t576 * t202 * t756 + t248 * (t203 * t206 * t765 + t760 * t759);
t777 = -t465 * t596 + t598;
t779 = t470 * t777 + t604;
t782 = t301 * t779;
t783 = t270 * t279;
t788 = -koppelP(5,3) * t465 - koppelP(5,1) * t694 + koppelP(5,2) * t696;
t799 = -t465 * t628 + t630;
t801 = t470 * t799 + t636;
t804 = t365 * t801;
t805 = t334 * t343;
t810 = -koppelP(6,3) * t465 - koppelP(6,1) * t694 + koppelP(6,2) * t696;
t822 = mrSges(4,1) * t462 - t663;
t833 = g(1) * mrSges(4,1);
t836 = -t257 * m(3) * (t203 * t759 - t215 * t765) + t608 * t266 * t779 + t312 * (t267 * t270 * t788 + t783 * t782) - t321 * m(3) * (t267 * t782 - t279 * t788) + t640 * t330 * t801 + t376 * (t331 * t334 * t810 + t805 * t804) - t385 * m(3) * (t331 * t804 - t343 * t810) - t470 * g(1) * (-t822 * t465 + mrSges(4,3) * t467) + t467 * t822 * g(3) - g(1) * mrSges(4,2) * t472 * t462 - t833 * t460 * t472 + mrSges(4,3) * t671;
t839 = t472 * t686 - t486;
t846 = t44 * t839;
t859 = t472 * t710 - t518;
t866 = t109 * t859;
t879 = t472 * t732 - t550;
t886 = t173 * t879;
t899 = t472 * t754 - t582;
t906 = t237 * t899;
t912 = t39 * (-t479 * t44 * t699 + t479 * t7 * t839) + t56 * (-t483 * t11 * t699 + t692 * t846) - t65 * m(3) * (t7 * t8 * t699 + t8 * t846) + t104 * (-t511 * t109 * t721 + t511 * t74 * t859) + t120 * (-t515 * t78 * t721 + t716 * t866) - t129 * m(3) * (t74 * t75 * t721 + t75 * t866) + t168 * (t543 * t138 * t879 - t543 * t173 * t743) + t184 * (-t547 * t142 * t743 + t738 * t886) - t193 * m(3) * (t138 * t139 * t743 + t139 * t886) + t232 * (t575 * t202 * t899 - t575 * t237 * t765) + t248 * (-t579 * t206 * t765 + t760 * t906);
t920 = t472 * t777 - t614;
t927 = t301 * t920;
t940 = t472 * t799 - t646;
t947 = t365 * t940;
t966 = t465 * t472;
t978 = -t257 * m(3) * (t202 * t203 * t765 + t203 * t906) + t296 * (t607 * t266 * t920 - t607 * t301 * t788) + t312 * (-t611 * t270 * t788 + t783 * t927) - t321 * m(3) * (t266 * t267 * t788 + t267 * t927) + t360 * (t639 * t330 * t940 - t639 * t365 * t810) + t376 * (-t643 * t334 * t810 + t805 * t947) - t385 * m(3) * (t330 * t331 * t810 + t331 * t947) - t467 * (mrSges(4,3) * t472 * g(1) + t462 * mrSges(4,1) * g(2) - mrSges(4,2) * t678) + t462 * (t470 * mrSges(4,2) + mrSges(4,1) * t966) * g(1) - mrSges(4,2) * g(1) * t460 * t966 + t833 * t460 * t470 - mrSges(4,3) * t658;
unknown(1,1) = t389;
unknown(2,1) = t433;
unknown(3,1) = t459;
unknown(4,1) = t684;
unknown(5,1) = t770 + t836;
unknown(6,1) = t912 + t978;
taugX  = unknown;
