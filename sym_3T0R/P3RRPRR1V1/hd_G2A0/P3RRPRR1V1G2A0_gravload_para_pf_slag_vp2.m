% Calculate Gravitation load for parallel robot
% P3RRPRR1V1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2022-11-04 17:08
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR1V1G2A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-04 17:08:05
% EndTime: 2022-11-04 17:08:05
% DurationCPUTime: 0.33s
% Computational Cost: add. (426->89), mult. (507->159), div. (36->7), fcn. (342->18), ass. (0->75)
t565 = legFrame(3,2);
t548 = sin(t565);
t551 = cos(t565);
t538 = t548 * g(1) + t551 * g(2);
t547 = pkin(1) * m(3) + mrSges(2,1) + mrSges(3,1);
t574 = cos(qJ(2,3));
t557 = 0.1e1 / t574;
t564 = mrSges(2,2) + mrSges(3,2);
t568 = sin(qJ(2,3));
t541 = t551 * g(1) - t548 * g(2);
t569 = sin(qJ(1,3));
t575 = cos(qJ(1,3));
t589 = g(3) * t575 + t541 * t569;
t611 = ((t538 * t564 + t589 * t547) * t568 - (t538 * t547 - t589 * t564) * t574) * t557;
t566 = legFrame(2,2);
t549 = sin(t566);
t552 = cos(t566);
t539 = t549 * g(1) + t552 * g(2);
t576 = cos(qJ(2,2));
t558 = 0.1e1 / t576;
t570 = sin(qJ(2,2));
t542 = t552 * g(1) - t549 * g(2);
t571 = sin(qJ(1,2));
t577 = cos(qJ(1,2));
t588 = g(3) * t577 + t542 * t571;
t610 = ((t539 * t564 + t588 * t547) * t570 - (t539 * t547 - t588 * t564) * t576) * t558;
t567 = legFrame(1,2);
t550 = sin(t567);
t553 = cos(t567);
t540 = t550 * g(1) + t553 * g(2);
t578 = cos(qJ(2,1));
t559 = 0.1e1 / t578;
t572 = sin(qJ(2,1));
t543 = t553 * g(1) - t550 * g(2);
t573 = sin(qJ(1,1));
t579 = cos(qJ(1,1));
t587 = g(3) * t579 + t543 * t573;
t609 = ((t540 * t564 + t587 * t547) * t572 - (t540 * t547 - t587 * t564) * t578) * t559;
t561 = pkin(3) + qJ(3,3);
t554 = 0.1e1 / t561;
t608 = (t569 * g(3) - t575 * t541) * t554;
t562 = pkin(3) + qJ(3,2);
t555 = 0.1e1 / t562;
t607 = (t571 * g(3) - t577 * t542) * t555;
t563 = pkin(3) + qJ(3,1);
t556 = 0.1e1 / t563;
t606 = (t573 * g(3) - t579 * t543) * t556;
t593 = mrSges(2,3) + mrSges(3,3) - mrSges(1,2);
t544 = m(3) * qJ(3,3) + t593;
t586 = -t547 * t574 + t564 * t568 - mrSges(1,1);
t605 = t554 * ((-t569 * t544 + t586 * t575) * t541 + (-t544 * t575 - t569 * t586) * g(3));
t545 = m(3) * qJ(3,2) + t593;
t585 = -t547 * t576 + t564 * t570 - mrSges(1,1);
t604 = t555 * ((-t571 * t545 + t585 * t577) * t542 + (-t545 * t577 - t571 * t585) * g(3));
t546 = m(3) * qJ(3,1) + t593;
t584 = -t547 * t578 + t564 * t572 - mrSges(1,1);
t603 = t556 * ((-t573 * t546 + t584 * t579) * t543 + (-t546 * t579 - t573 * t584) * g(3));
t580 = pkin(2) + pkin(1);
t602 = t568 * t580;
t601 = t569 * t574;
t600 = t570 * t580;
t599 = t571 * t576;
t598 = t572 * t580;
t597 = t573 * t578;
t596 = t574 * t580;
t595 = t576 * t580;
t594 = t578 * t580;
t592 = t557 * t605;
t591 = t558 * t604;
t590 = t559 * t603;
t583 = -t561 * t575 + t569 * t596;
t582 = -t562 * t577 + t571 * t595;
t581 = -t563 * t579 + t573 * t594;
t560 = 0.1e1 / t580;
t1 = [(t550 * t572 + t553 * t597) * t590 + (t549 * t570 + t552 * t599) * t591 + (t548 * t568 + t551 * t601) * t592 - g(1) * m(4) + (t548 * t611 + t549 * t610 + t550 * t609) * t560 + (-(t550 * t598 + t581 * t553) * t606 - (t549 * t600 + t582 * t552) * t607 - (t548 * t602 + t583 * t551) * t608) * m(3); (-t550 * t597 + t572 * t553) * t590 + (-t549 * t599 + t570 * t552) * t591 + (-t548 * t601 + t568 * t551) * t592 - g(2) * m(4) + (t551 * t611 + t552 * t610 + t553 * t609) * t560 + (-(-t581 * t550 + t553 * t598) * t606 - (-t582 * t549 + t552 * t600) * t607 - (-t583 * t548 + t551 * t602) * t608) * m(3); t575 * t605 + t577 * t604 + t579 * t603 - g(3) * m(4) + (-(t573 * t563 + t579 * t594) * t606 - (t571 * t562 + t577 * t595) * t607 - (t569 * t561 + t575 * t596) * t608) * m(3);];
taugX  = t1;
