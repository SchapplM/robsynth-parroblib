% Calculate Gravitation load for parallel robot
% P3PRRR1G3P3A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2020-03-09 21:07
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRR1G3P3A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:06:51
% EndTime: 2020-03-09 21:06:52
% DurationCPUTime: 0.29s
% Computational Cost: add. (474->68), mult. (399->106), div. (45->5), fcn. (291->18), ass. (0->60)
t582 = legFrame(3,2);
t571 = sin(t582);
t574 = cos(t582);
t552 = t574 * g(1) - t571 * g(2);
t578 = pkin(7) + qJ(2,3);
t568 = qJ(3,3) + t578;
t555 = sin(t568);
t558 = cos(t568);
t597 = mrSges(3,2) * g(3);
t598 = mrSges(3,1) * g(3);
t531 = (mrSges(3,2) * t552 + t598) * t558 + (mrSges(3,1) * t552 - t597) * t555;
t562 = sin(t578);
t565 = cos(t578);
t534 = 0.1e1 / (t555 * t565 - t562 * t558);
t604 = t531 * t534;
t583 = legFrame(2,2);
t572 = sin(t583);
t575 = cos(t583);
t553 = t575 * g(1) - t572 * g(2);
t579 = pkin(7) + qJ(2,2);
t569 = qJ(3,2) + t579;
t556 = sin(t569);
t559 = cos(t569);
t532 = (mrSges(3,2) * t553 + t598) * t559 + (mrSges(3,1) * t553 - t597) * t556;
t563 = sin(t579);
t566 = cos(t579);
t535 = 0.1e1 / (t556 * t566 - t563 * t559);
t603 = t532 * t535;
t584 = legFrame(1,2);
t573 = sin(t584);
t576 = cos(t584);
t554 = t576 * g(1) - t573 * g(2);
t580 = pkin(7) + qJ(2,1);
t570 = qJ(3,1) + t580;
t557 = sin(t570);
t560 = cos(t570);
t533 = (mrSges(3,2) * t554 + t598) * t560 + (mrSges(3,1) * t554 - t597) * t557;
t564 = sin(t580);
t567 = cos(t580);
t536 = 0.1e1 / (t557 * t567 - t564 * t560);
t602 = t533 * t536;
t601 = (pkin(2) * t565 + pkin(3) * t558) * t604;
t600 = (pkin(2) * t566 + pkin(3) * t559) * t603;
t599 = (pkin(2) * t567 + pkin(3) * t560) * t602;
t596 = g(3) * mrSges(2,2);
t577 = m(3) * pkin(2) + mrSges(2,1);
t561 = g(3) * t577;
t592 = t534 * ((mrSges(2,2) * t552 + t561) * t565 + t562 * (t577 * t552 - t596) + t531);
t591 = t535 * ((mrSges(2,2) * t553 + t561) * t566 + t563 * (t577 * t553 - t596) + t532);
t590 = t536 * ((mrSges(2,2) * t554 + t561) * t567 + t564 * (t577 * t554 - t596) + t533);
t589 = t558 * t592;
t588 = t559 * t591;
t587 = t560 * t590;
t586 = 0.1e1 / pkin(2);
t585 = 0.1e1 / pkin(3);
t581 = m(1) + m(2) + m(3);
t551 = t573 * g(1) + t576 * g(2);
t550 = t572 * g(1) + t575 * g(2);
t549 = t571 * g(1) + t574 * g(2);
t1 = [-g(1) * m(4) + (-t571 * t549 - t572 * t550 - t573 * t551) * t581 + (t574 * t589 + t575 * t588 + t576 * t587 + (-t574 * t601 - t575 * t600 - t576 * t599) * t585) * t586; -g(2) * m(4) + (-t549 * t574 - t550 * t575 - t551 * t576) * t581 + (-t571 * t589 - t572 * t588 - t573 * t587 + (t571 * t601 + t572 * t600 + t573 * t599) * t585) * t586; -g(3) * m(4) + (-t555 * t592 - t556 * t591 - t557 * t590 + ((pkin(2) * t564 + pkin(3) * t557) * t602 + (pkin(2) * t563 + pkin(3) * t556) * t603 + (pkin(2) * t562 + pkin(3) * t555) * t604) * t585) * t586;];
taugX  = t1;
