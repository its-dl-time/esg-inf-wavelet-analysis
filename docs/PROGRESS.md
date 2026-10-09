# Nhật ký công việc đến 2026-10-09

| Giai đoạn | Trạng thái | Có thể kiểm tra trong repo |
|---|---|---|
| Dữ liệu EPU/CPI, 18 nước, mùa vụ STL log | Đã chạy ở pipeline gốc | `data/`, `R/prepare_data.R`, `data/processed/seasonality_scores.csv` |
| WTC ban đầu, 149 công việc, kiểm định nhóm chu kỳ, chẩn đoán | Đã chạy với smoother 5/0,25 | `robustness/legacy/`, `results/legacy/` |
| Pilot ba nước, ba null, smoother 5/0,25 và 9/1 | Đã được tóm tắt trong cuộc trao đổi trước; các ma trận/checkpoint của lượt đó không còn trong workspace khi đóng gói | `docs/RESULTS.md` ghi lại số liệu đã báo, cần chạy lại mã mới để xác nhận |
| Lõi WTC 9/1 cho 18 nước | Mã và dữ liệu đã chuẩn bị; chưa chạy toàn bộ 18 nước | `R/core_wtc.R`, đầu ra được tạo trong `results/core/` |
| UCT và Brent, partial/multiple | Mã mở rộng; đầu vào ancillary tải từ nguồn khi chạy; không có p-value điều kiện hiệu chuẩn | `extensions/uct/` |

Môi trường đóng gói lần này không có Rscript, vì vậy chưa thể thực thi lại mã R mới tại đây. Bảng cũ và tóm tắt pilot không được dán nhãn là kết quả của lõi mới. Kiểm tra tĩnh và tính toàn vẹn dữ liệu được ghi trong manifest.
