# Kết quả đã có và giới hạn diễn giải

## 18 nước: bản ước lượng cũ

`results/legacy/main_band_summary.csv` có 54 tổ hợp nước×nhóm chu kỳ từ WaveletComp với smoother thời gian 5, scale 0,25 và null AR(1). Kiểm định mean nhóm chu kỳ nhỏ nhất là Canada 8–<32 tháng, p = 0,046953; **0/54** qua Holm ở 5%. Các kiểm định BH/BY ở cấp từng ô WTC của lượt gốc cũng không giữ ô có ý nghĩa. Bản đồ cũ nhạy với bộ ước lượng: trung vị tương quan bản đồ WaveletComp và biwavelet khoảng 0,349 trên vùng chung, theo báo cáo gốc. Đây không phải kết quả của smoother 9/1.

## Pilot ba nước đã báo trước khi đóng gói

Ở smoother 9/1, null AR(1), Fourier và AR(1) phương sai cục bộ với 250 map hiệu chuẩn + 999 map kiểm định: **0/9** nước×nhóm chu kỳ qua Holm riêng cho cả band mean, max cluster và max episode dưới từng null. Các cluster dài hạn có p thô xấp xỉ Mỹ 0,044/0,099/0,030 và Trung Quốc 0,044/0,043/0,041 (AR/Fourier/local-volatility); Holm tương ứng Mỹ 0,352/0,595/0,270 và Trung Quốc 0,352/0,387/0,328. Canada cụm ngắn hạn có p thô 0,036/0,046/0,132. Những số này là bản tóm tắt từ lượt trước, **không có checkpoint pilot để kiểm lại tại thời điểm tạo repo**.

Chọn 9/1 chủ yếu vì bản đồ gần biwavelet hơn trong pilot (tương quan Mỹ 0,759; Trung Quốc 0,847; Canada 0,796) và ma trận conditional tại dải dài ổn định hơn. Tương quan là đo độ nhạy, không tự nó chứng minh hai bộ ước lượng có cùng suy luận. Cần chạy `R/core_wtc.R` với ba nước trước để tái lập, sau đó 18 nước với cùng phiên bản code.

UCT chính thức trong pilot trước kết thúc 2024-02; cửa sổ khớp khi ấy là Mỹ 1993-02–2024-02 (373 tháng) và Trung Quốc 1997-02–2024-02 (325 tháng). Partial đôi khi tăng WTC; multiple gần 1 khi thêm nhiều dự báo, đặc biệt khi ma trận gần suy biến. Đây là mô tả, không suy luận p về hiệu ứng riêng EPU. Extension mới dùng biwavelet và biến đổi UCT được khai báo rõ; kết quả không được trộn số học với pilot cũ nếu dữ liệu/định nghĩa khác.

## Quy tắc diễn giải khi chạy toàn bộ

Đọc `results/core/summary.csv` theo từng `statistic`. `p` hiệu chỉnh sự tìm kiếm *bên trong* từng bản đồ nước×band qua phân phối max; `q_BH` xử lý family country×band của chính thống kê đó; `q_BY` là đối chiếu khi không giả định phụ thuộc dương; `p_Holm` là đối chiếu kiểm soát FWER. Không chọn một trong ba thống kê hay null có p nhỏ nhất sau khi xem kết quả rồi gọi là phép thử duy nhất. Nếu ba thống kê đều được coi là phát hiện độc lập chính thức, family còn phải mở rộng hoặc đặt trước thứ tự kiểm định. Cú sốc nhân quả cần định nghĩa biến cố và mô hình nhận dạng riêng; góc pha không đủ.
